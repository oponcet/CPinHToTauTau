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
import matplotlib.pyplot as plt


logger = law.logger.get_logger(__name__)


def add_hist_hooks(config: od.Config) -> None:
    """
    Add histogram hooks to a configuration.
    """
    print("add_fake_factors")
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

    def qcd_estimation(task, hists):
        print(f"qcd_estimation")
        # If `hists` is empty, return it immediately.
        
        if not hists:
            print("no hists")
            return hists

        # Get the QCD process from the config. If it doesn't exist, return the histograms unchanged.
        qcd_proc = config.get_process("qcd", default=None)
        if not qcd_proc:
            return hists

        # extract all unique category ids and verify that the axis order is exactly
        # "category -> shift -> variable" which is needed to insert values at the end
        CAT_AXIS, SHIFT_AXIS, VAR_AXIS = range(3)

        # Collect all unique category IDs from the histograms and validate axis order
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

        # Create QCD groups (os_iso, os_noniso, ss_iso, ss_noniso) based on category tags
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        for cat_id in category_ids:
            # Get the category instance using its ID
            cat_inst = config.get_category(cat_id)
            # Populate the QCD groups by checking specific category tags
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_iso = cat_inst
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].os_noniso = cat_inst
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_iso = cat_inst
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups[cat_inst.x.qcd_group].ss_noniso = cat_inst

        # Filter out QCD groups that contain all four required regions (os_iso, os_noniso, ss_iso, ss_noniso)
        complete_groups = [name for name, cats in qcd_groups.items() if len(cats) == 4]

        # nothing to do if there are no complete groups
        if not complete_groups:
            return hists

        # Sum up the Monte Carlo (MC) histograms and data histograms. Return early if either is empty.
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        # start by copying the mc hist and reset it, then fill it at specific category slices
        hists[qcd_proc] = qcd_hist = mc_hist.copy().reset()
        for group_name in complete_groups:
            group = qcd_groups[group_name]

            # get the corresponding histograms and convert them to number objects,
            # each one storing an array of values with uncertainties
            # shapes: (SHIFT, VAR)
            get_hist = lambda h, region_name: h[{"category": hist.loc(group[region_name].id)}]
            os_noniso_mc = hist_to_num(get_hist(mc_hist, "tautau_antiIso"), "tautau_antiIso_mc") # is_noniso
            ss_noniso_mc = hist_to_num(get_hist(mc_hist, "FFDRantiIso_tautau"), "FFDRantiIso_tautau_mc") # ss_noniso
            ss_iso_mc = hist_to_num(get_hist(mc_hist, "FFDRIso_tautau"), "FFDRIso_tautau_mc") # ss_iso
            os_noniso_data = hist_to_num(get_hist(data_hist, "tautau_antiIso"), "tautau_antiIso_data")
            ss_noniso_data = hist_to_num(get_hist(data_hist, "FFDRantiIso_tautau"), "FFDRantiIso_tautau_data")
            ss_iso_data = hist_to_num(get_hist(data_hist, "FFDRantiIso_tautau_iso"), "FFDRantiIso_tautau_data")

            # estimate qcd shapes in the three sideband regions
            # shapes: (SHIFT, VAR)
            os_noniso_qcd = os_noniso_data - os_noniso_mc
            ss_iso_qcd = ss_iso_data - ss_iso_mc
            ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

            # get integrals in ss regions for the transfer factor
            # shapes: (SHIFT,)
            int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
            int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)

            # complain about negative integrals
            int_ss_iso_neg = int_ss_iso <= 0
            int_ss_noniso_neg = int_ss_noniso <= 0
            if int_ss_iso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_iso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )
            if int_ss_noniso_neg.any():
                shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_ss_noniso_neg)[0]))
                shifts = list(map(config.get_shift, shift_ids))
                logger.warning(
                    f"negative QCD integral in ss_noniso region for group {group_name} and shifts: "
                    f"{', '.join(map(str, shifts))}",
                )

            # ABCD method
            # shape: (SHIFT, VAR)
            os_iso_qcd = os_noniso_qcd * ((int_ss_iso / int_ss_noniso)[:, None])

            # combine uncertainties and store values in bare arrays
            os_iso_qcd_values = os_iso_qcd()
            os_iso_qcd_variances = os_iso_qcd(sn.UP, sn.ALL, unc=True)**2

            # define uncertainties
            unc_data = os_iso_qcd(sn.UP, ["os_noniso_data", "ss_iso_data", "ss_noniso_data"], unc=True)
            unc_mc = os_iso_qcd(sn.UP, ["os_noniso_mc", "ss_iso_mc", "ss_noniso_mc"], unc=True)
            unc_data_rel = abs(unc_data / os_iso_qcd_values)
            unc_mc_rel = abs(unc_mc / os_iso_qcd_values)

            # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
            keep_variance_mask = (
                np.isfinite(unc_mc_rel) &
                (unc_mc_rel > unc_data_rel) &
                (unc_mc_rel > 0.15)
            )
            os_iso_qcd_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
            os_iso_qcd_variances[~keep_variance_mask] = 0

            # retro-actively set values to zero for shifts that had negative integrals
            neg_int_mask = int_ss_iso_neg | int_ss_noniso_neg
            os_iso_qcd_values[neg_int_mask] = 1e-5
            os_iso_qcd_variances[neg_int_mask] = 0

            # residual zero filling
            zero_mask = os_iso_qcd_values <= 0
            os_iso_qcd_values[zero_mask] = 1e-5
            os_iso_qcd_variances[zero_mask] = 0

            # insert values into the qcd histogram
            cat_axis = qcd_hist.axes["category"]
            for cat_index in range(cat_axis.size):
                if cat_axis.value(cat_index) == group.os_iso.id:
                    qcd_hist.view().value[cat_index, ...] = os_iso_qcd_values
                    qcd_hist.view().variance[cat_index, ...] = os_iso_qcd_variances
                    break
            else:
                raise RuntimeError(
                    f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                    f"for category {group.os_iso}",
                )

        return hists


    def fake_factor(task, hists):
        print(f"fake_factor")

        if not hists:
            print("no hists")
            return hists

        # get dummy processes
        factor_bin = config.get_process("qcd", default=None)
        print(f"factor_bin: {factor_bin}")

        if not factor_bin:
            return hists

        # factor_int = config.get_process("dy_lep_m50", default=None)
        # print(f"factor_int: {factor_int}")
        # if not factor_int:
        #     return hists

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

        print("category_ids: ", category_ids)
        # create qcd groups of Signal region 
        qcd_groups: dict[str, dict[str, od.Category]] = defaultdict(DotDict)
        print("qcd_groups: ", qcd_groups)
        for cat_id in category_ids:
            print("cat_id: ", cat_id)
            cat_inst = config.get_category(cat_id)
            print("cat_inst: ", cat_inst)
            if cat_inst.has_tag({"os", "iso"}, mode=all):
                qcd_groups["tautau"] = cat_inst  # Signal region
                print("category id: ", cat_inst.id , "is os_iso")
            elif cat_inst.has_tag({"os", "noniso"}, mode=all):
                qcd_groups["tautau_antiIso"] = cat_inst  # Application region
                print("category id: ", cat_inst.id , "is os_iso")
            elif cat_inst.has_tag({"ss", "iso"}, mode=all):
                qcd_groups["FFDRIso_tautau"] = cat_inst  # Derivative region
                print("category id: ", cat_inst.id , "is os_iso")
            elif cat_inst.has_tag({"ss", "noniso"}, mode=all):
                qcd_groups["FFDRantiIso_tautau"] = cat_inst  # Derivative region
                print("category id: ", cat_inst.id , "is os_iso")

        # get complete qcd groups
        complete_groups = [name for name, cats in qcd_groups.items() if len(qcd_groups) == 4]
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
        print("mc_hist: ", mc_hist)
        print("data_hist: ", data_hist)

        # start by copying the mc hist and reset it, then fill it at specific category slices
        # hists = {}
        hists[factor_bin] = factor_hist = mc_hist.copy().reset()
        # hists[factor_int] = factor_hist_int = mc_hist.copy().reset()
        


        # get the corresponding histograms and convert them to number objects,
        # each one storing an array of values with uncertainties
        # shapes: (SHIFT, VAR)
        # from IPython import embed; embed()
        get_hist = lambda h, region_name: h[{"category": hist.loc(qcd_groups[region_name].id)}]
        ss_noniso_mc = hist_to_num(get_hist(mc_hist, "FFDRantiIso_tautau"), "FFDRantiIso_tautau_mc") # ss_noniso
        ss_iso_mc = hist_to_num(get_hist(mc_hist, "FFDRIso_tautau"), "FFDRIso_tautau_mc") # ss_iso
        ss_noniso_data = hist_to_num(get_hist(data_hist, "FFDRantiIso_tautau"), "FFDRantiIso_tautau_data") # ss_noniso
        ss_iso_data = hist_to_num(get_hist(data_hist, "FFDRIso_tautau"), "FFDRIso_tautau_data") # ss_iso

        # take the difference between data and MC in the control regions
        ss_iso_qcd = ss_iso_data - ss_iso_mc
        ss_noniso_qcd = ss_noniso_data - ss_noniso_mc

        # # calculate the pt-independent fake factor
        # int_ss_iso = integrate_num(ss_iso_qcd, axis=1)
        # int_ss_noniso = integrate_num(ss_noniso_qcd, axis=1)
        # fake_factor_int = (int_ss_iso / int_ss_noniso)[0, None]

        # calculate the pt-dependent fake factor
        fake_factor = (ss_iso_qcd / ss_noniso_qcd)[:, None]
        fake_factor_values = np.squeeze(np.nan_to_num(fake_factor()), axis=0)
        fake_factor_variances = fake_factor(sn.UP, sn.ALL, unc=True)**2

        # print(f"fake_factor_values: {fake_factor_values}")
        # print(f"fake_factor_variances: {fake_factor_variances}")

        # # change shape of fake_factor_int for plotting
        # fake_factor_int_values = fake_factor_values.copy()
        # fake_factor_int_values.fill(fake_factor_int()[0])

        # srint(f"fake_factor_int_values: {fake_factor_int_values}")

        # plot the fake factor
        # Plot the fake factor with error bars
        # pt_bins = factor_hist.axes["tau_1_pt"].edges  # Example of how to extract bin edges (replace 'pt' with actual axis name)
        # pt_centers = (pt_bins[:-1] + pt_bins[1:]) / 2  # Get the bin centers for plotting

        # fake_factor_uncertainties = np.sqrt(fake_factor_variances)  # Calculate uncertainties

        # # Plot using Matplotlib
        # plt.figure(figsize=(8, 6))
        # plt.errorbar(pt_centers, fake_factor_values, yerr=fake_factor_variances, fmt='o', label='Fake Factor', color='b')
        # plt.xlabel('pT [GeV]')  # Adjust x-axis label as necessary
        # plt.ylabel('Fake Factor')
        # plt.title(f'Fake Factor vs pT for group {group_name}')
        # plt.grid(True)
        # plt.legend()
        # plt.saveAs("FakeFactors.png")  # save


        # insert values into the qcd histogram
        cat_axis = factor_hist.axes["category"]
        print(f"cat_axis: {cat_axis}")
        for cat_index in range(cat_axis.size):
            print(f"cat_axis.value(cat_index): {cat_axis.value(cat_index)}")
            if cat_axis.value(cat_index) == qcd_groups["tautau"].id:
                print(f"cat_axis.value(cat_index): {cat_axis.value(cat_index)}")
                # print(f"fake factor for category tautau: {fake_factor_int()}")
                factor_hist.view().value[cat_index, ...] = fake_factor_values
                factor_hist.view().variance[cat_index, ...] = fake_factor_variances
                # factor_hist_int.view().value[cat_index, ...] = fake_factor_int_values
                print(f"Fake factors values = {fake_factor_values}")
                break
        else:
            raise RuntimeError(
                f"could not find index of bin on 'category' axis of qcd histogram {factor_hist} "
                f"for category tautau ",
            )

        print("hist: ", hists)
        return hists

    
    def example_hook(task, hists):
    # create a new "QCD" process, if not already done in the config itself
        print("example_hook")

        import order as od
        qcd = od.Process("qcd", id="+", label="QCD", color1=(244, 93, 66))

        from IPython import embed; embed()

        # Sum up the Monte Carlo (MC) histograms and data histograms. Return early if either is empty.
        mc_hists = [h for p, h in hists.items() if p.is_mc and not p.has_tag("signal")]
        data_hists = [h for p, h in hists.items() if p.is_data]
        if not mc_hists or not data_hists:
            return hists
        mc_hist = sum(mc_hists[1:], mc_hists[0].copy())
        data_hist = sum(data_hists[1:], data_hists[0].copy())

        hists[qcd] = mc_hist - data_hist

        mc_hist_num = hist_to_num(mc_hist, "mc")
        data_hist_num = hist_to_num(data_hist, "data")

        qcd_num = mc_hist_num - data_hist_num

        # complain about negative integrals
        int_qcd_num_neg = qcd_num <= 0
        if int_qcd_num_neg.any():
            shift_ids = list(map(mc_hist.axes["shift"].value, np.where(int_qcd_num_neg)[0]))
            shifts = list(map(config.get_shift, shift_ids))
            logger.warning(
                f"negative QCD integral in ss_iso region for group {group_name} and shifts: "
                f"{', '.join(map(str, shifts))}",
            )
        

        # combine uncertainties and store values in bare arrays
        qcd_num_values = qcd_num()
        qcd_num_variances = qcd_num(sn.UP, sn.ALL, unc=True)**2

        # define uncertainties
        unc_data = qcd_num(sn.UP, ["data"], unc=True)
        unc_mc = qcd_num(sn.UP, ["mc"], unc=True)
        unc_data_rel = abs(unc_data / qcd_num_values)
        unc_mc_rel = abs(unc_mc / qcd_num_values)
                

        # only keep the MC uncertainty if it is larger than the data uncertainty and larger than 15%
        keep_variance_mask = (
            np.isfinite(unc_mc_rel) &
            (unc_mc_rel > unc_data_rel) &
            (unc_mc_rel > 0.15)
        )
        qcd_num_variances[keep_variance_mask] = unc_mc[keep_variance_mask]**2
        qcd_num_variances[~keep_variance_mask] = 0

        # retro-actively set values to zero for shifts that had negative integrals
        neg_int_mask = int_qcd_num_neg
        qcd_num_values[neg_int_mask] = 1e-5
        qcd_num_variances[neg_int_mask] = 0

        # residual zero filling
        zero_mask = qcd_num_values <= 0
        qcd_num_values[zero_mask] = 1e-5
        qcd_num_variances[zero_mask] = 0

        # insert values into the qcd histogram
        cat_axis = qcd_hist.axes["category"]
        for cat_index in range(cat_axis.size):
            if cat_axis.value(cat_index) == group.os_iso.id:
                qcd_hist.view().value[cat_index, ...] = qcd_num_values
                qcd_hist.view().variance[cat_index, ...] = qcd_num_variances
                break
        else:
            raise RuntimeError(
                f"could not find index of bin on 'category' axis of qcd histogram {qcd_hist} "
                f"for category {group.os_iso}",
            )

        return hists
    
    config.x.hist_hooks = {
        "example_hook": example_hook,
        "qcd_estimation": qcd_estimation,
        "fake_factor": fake_factor,
    }