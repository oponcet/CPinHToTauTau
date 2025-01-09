import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

from httcp.util import get_trigger_id_map

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
warn   = maybe_import("warnings")

# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


# ------------------------------------------------- #
# Assign MC weight [gen Weight / LHEWeight]
# ------------------------------------------------- #
@producer(
    uses={"genWeight", optional("LHEWeight.originalXWGTUP")},
    produces={"mc_weight"},
    mc_only=True,
)
def scale_mc_weight(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    """
    Reads the genWeight and LHEWeight columns and makes a decision about which one to save. This
    should have been configured centrally [1] and stored in genWeight, but there are some samples
    where this failed.

    Strategy:

      1. Use LHEWeight.originalXWGTUP when it exists and genWeight is always 1.
      2. In all other cases, use genWeight.

    [1] https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD?rev=99#Weigths
    """
    # determine the mc_weight
    mc_weight = np.sign(events.genWeight)
    if has_ak_column(events, "LHEWeight.originalXWGTUP") and ak.all(events.genWeight == 1.0):
        mc_weight = np.sign(events.LHEWeight.originalXWGTUP)

    # store the column
    events = set_ak_column(events, "mc_weight", mc_weight, value_type=np.float32)

    return events


# ------------------------------------------------- #
# Calculate Zpt weights
# ------------------------------------------------- #

@producer(
    uses={
        "GenZ.pt", "GenZ.mass",
    },
    produces={
        "zpt_reweight",
    },
    mc_only=True,
)
def zpt_reweight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # if no gen particle found, all fields of GenZ will be zero
    # Those events will have nominal zpt rewt = 1.0
    # the same will be applied also for the evnts outside the range

    is_outside_range = (
        ((events.GenZ.pt == 0.0) & (events.GenZ.mass == 0.0))
        | ((events.GenZ.pt >= 600.0) | (events.GenZ.mass >= 1000.0))
    )

    # for safety
    zm  = ak.where(events.GenZ.mass > 1000.0, 999.99, events.GenZ.mass)
    zpt = ak.where(events.GenZ.pt > 600.0, 599.99, events.GenZ.pt)

    sf_nom = ak.where(is_outside_range,
                      1.0,
                      self.zpt_corrector.evaluate(zm,zpt))

    events = set_ak_column(events, "zpt_reweight", sf_nom, value_type=np.float32)
    
    return events


@zpt_reweight.requires
def zpt_reweight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@zpt_reweight.setup
def zpt_reweight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate
    
    correction_set = correctionlib.CorrectionSet.from_string(
        bundle.files.zpt_rewt_sf.load(formatter="gzip").decode("utf-8"),
    )
    self.zpt_corrector    = correction_set["zptreweight"]





# ------------------------------------------------- #
# Calculate FF weights (Dummy)
# ------------------------------------------------- #

@producer(
    uses={
        "channel_id",
        "category_ids",
        "hcand.pt", "hcand.mass", "PuppiMET.pt", "PuppiMET.phi",
    },
    produces={
        "ff_weight",
        "closure_weight"
    },
    mc_only=False,
)
def ff_weight(
        self: Producer,
        events: ak.Array,
        **kwargs,
) :
    # Leading candidate
    from IPython import embed; embed()
    hcand1 = events.hcand[:,0:1] # Leading candidate hcand1 = events.hcand[:,0] 


    # If events.is_os, events.is_real_1, events.is_iso_2 are True and events.is_iso_1 is False, then is_C_category is True
    is_C_category = (
    events.is_os &        # All of these conditions must be True
    # events.is_real_1 &
    events.is_iso_2 &
    ~events.is_iso_1      # This must be False (negate with ~)
    )

    is_C_category = flat_np_view(is_C_category)

    # Get the pt of the leading candidate
    pt1 = flat_np_view(hcand1.pt[:,None])

    # Get met_var_qcd_h1
    met = ak.with_name(events.PuppiMET, "PtEtaPhiMLorentzVector")
    dphi_met_h1 = met.delta_phi(hcand1)
    met_var_qcd_h1 = met.pt * np.cos(dphi_met_h1)/hcand1.pt
    met_var_qcd_h1 = flat_np_view(met_var_qcd_h1[:,None])


    # Get the decay mode and jet multiplicity
    dms = ["a1dm11_1", "a1dm10_1", "a1dm2_1", "pi_1", "rho_1"]  # Decay modes
    njets = ["has_0j", "has_1j", "has_2j"]  # Jet multiplicity

    for dm in dms:
        for njet in njets:
            fake_factors_nom = self.ff_corrector.evaluate(
                pt1,
                dm,
                njet
            )
            closure_nom = self.closure_corrector.evaluate(
                met_var_qcd_h1,
                dm,
                njet
            )

    # Apply the fake factor only for the C category
    ff_nom = np.where(is_C_category, fake_factors_nom, 1.0)

    closure_nom = np.where(is_C_category, closure_nom, 1.0)


    # Add the column to the events
    events = set_ak_column(events, "ff_weight", ff_nom, value_type=np.float32)
    events = set_ak_column(events, "closure_weight", closure_nom, value_type=np.float32)
    
    return events


@ff_weight.requires
def ff_weight_requires(self: Producer, reqs: dict) -> None:
    if "external_files" in reqs:
        return
    
    from columnflow.tasks.external import BundleExternalFiles
    reqs["external_files"] = BundleExternalFiles.req(self.task)

@ff_weight.setup
def ff_weight_setup(
    self: Producer,
    reqs: dict,
    inputs: dict,
    reader_targets: InsertableDict,
) -> None:
    bundle = reqs["external_files"]
    import correctionlib
    correctionlib.highlevel.Correction.__call__ = correctionlib.highlevel.Correction.evaluate

    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.tautau_ff.path,
    ) 

    self.ff_corrector    = correction_set["fake_factors_fit"]

    correction_set = correctionlib.CorrectionSet.from_file(
        bundle.files.tautau_closure.path,
    ) 

    self.closure_corrector    = correction_set["closure_correction_fit"]
