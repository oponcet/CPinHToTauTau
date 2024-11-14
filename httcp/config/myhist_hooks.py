# coding: utf-8

"""
Histogram hooks.
"""

from __future__ import annotations

from collections import defaultdict

import law
import order as od
import scinum as sn

from columnflow.util import maybe_import, DotDict

np = maybe_import("numpy")
hist = maybe_import("hist")


logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    # helper to convert a histogram to a number object containing bin values and uncertainties
    # from variances stored in an array of values
    def hist_to_num(h: hist.Histogram, unc_name=str(sn.DEFAULT)) -> sn.Number:
        # Return an `sn.Number` object, where:
        # - The nominal values (i.e., bin contents) are `h.values()`.
        # - The uncertainties are stored in a dictionary where the key is `unc_name` and the value 
        #   is the square root of the variances (h.variances()**0.5).
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})
        return sn.Number(h.values(), {unc_name: h.variances()**0.5})

    # helper to integrate values stored in an array based number object
    # Define a helper function `integrate_num` that integrates the values of an `sn.Number` object `num`.
    # `axis` is an optional argument for specifying which axis to integrate along (if applicable).
    
    def integrate_num(num: sn.Number, axis=None) -> sn.Number:
        # Return a new `sn.Number` object, where:
        # - `nominal` represents the sum of the nominal values along the specified `axis`.
        # - The uncertainties are calculated for each `unc_name` in the `num.uncertainties` dictionary.
        #   The uncertainty values are computed by summing the squares of the upper and lower uncertainties 
        #   (stored as `unc_values_up` and `unc_values_down`), and taking the square root of the sum.
        
        return sn.Number(
            nominal=num.nominal.sum(axis=axis), # Sum the nominal values along the specified axis.
            uncertainties={
                unc_name: (
                    (unc_values_up**2).sum(axis=axis)**0.5,
                    (unc_values_down**2).sum(axis=axis)**0.5,
                )
                for unc_name, (unc_values_up, unc_values_down) in num.uncertainties.items()
            },
        )


    def fake_factor(task, hists):
        if not hists:
            return hists

        # get dummy processes
        factor_bin = config.get_process("qcd", default=None)
        if not factor_bin:
            return hists

        factor_int = config.get_process("dy", default=None)
        if not factor_int:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)
        category_ids = set()
        for proc, h in hists.items():
            # validate axes
            assert len(h.axes) == 3
            assert h.axes[CAT_AXIS].name == "category"
            assert h.axes[SHIFT_AXIS].name == "shift"
            # get the category axis
            cat_ax = h.axes["category"]
            for cat_index in range(cat_ax.size):
                category_ids.add(cat_ax.value(cat_index))

        # create qcd groups of Signal region 
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            cat_inst = config.get_category(cat_id)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # sum up mc and data histograms, stop early when empty
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists = {}
        hists[factor_bin] = factor_hist = mc_hist.copy().reset()
        hists[factor_int] = factor_hist_int = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]
            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "ss_noniso"), "ss_noniso_mc")
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "ss_iso"), "ss_iso_mc")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "ss_noniso"), "ss_noniso_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "ss_iso"), "ss_iso_data")

            # take the difference between data and MC in the control regions
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # calculate the pt-independent fake factor
            int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)
            fake_factor_int = (int_ss_iso / int_ss_noniso)[0, None]

            # calculate the pt-dependent fake factor
            fake_factor = (ss_iso_qcd / ss_noniso_qcd)[:, None]
            fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0)
            fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2


            # change shape of fake_factor_int for plotting
            fake_factor_int_values = fake_factor_values.copy()
            fake_factor_int_values.fill(fake_factor_int()[0])

            # insert values into the qcd histogram
            cat_axis = factor_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    factor_hist.view().value[cat_index, ...] = fake_factor_values
                    factor_hist.view().variance[cat_index, ...] = fake_factor_variances
                    factor_hist_int.view().value[cat_index, ...] = fake_factor_int_values
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {factor_hist} "
                    f"for category {group.os_iso}",
                )
        return hists

    config.x.hist_hooks = {
        "fake_factor": fake_factor,
    }