import functools
from columnflow.production import Producer, producer
from columnflow.util import maybe_import, safe_div, InsertableDict
from columnflow.columnar_util import set_ak_column,remove_ak_column, has_ak_column, EMPTY_FLOAT, Route, flat_np_view, optional_column as optional
from columnflow.production.util import attach_coffea_behavior

ak     = maybe_import("awkward")
np     = maybe_import("numpy")
coffea = maybe_import("coffea")
# helper
set_ak_column_f32 = functools.partial(set_ak_column, value_type=np.float32)


@producer(
    uses={f"Tau.{var}" for var in [
                "decayMode", "genPartFlav"
                ] 
    } | {"process_id"},
    produces={
        "process_id"
    },
    mc_only=True,
)
def split_dy(self: Producer, events: ak.Array, **kwargs) -> ak.Array:
    

    tau_part_flav = {
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5
    }
    hcand_ele_mu_DM = events.hcand.decayMode[:,0]
    # hcand_Tau_idx = events.hcand.rawIdx[:,1:2]
    # mask_hcand_tau_idx = ak.flatten(hcand_Tau_idx)
    match = events.Tau.genPartFlav
    #### The ak.any() is a nasty fix that need to be removed, we have tautau that needs to be taken into account properly
    mu2tau_fakes_mask = ak.any(((match == tau_part_flav["prompt_mu"]) | (match == tau_part_flav["tau->mu"])),axis=1)
    e2tau_fakes_mask = ak.any(((match == tau_part_flav["prompt_e"]) | (match == tau_part_flav["tau->e"])),axis=1)
    genuine_tau_mask = ak.any((match == tau_part_flav["tau_had"]),axis=1)

    z_proc_id = self.dataset_inst.processes.values()[0].id
    ztoee_proc_id = 51101
    ztomm_proc_id = 51102
    ztott_proc_id = 51104
    ztoll_proc_id = 51103

    #from IPython import embed; embed()
    
    #process_id = np.array(events.process_id, dtype=np.int64)

    ztoeeormm = (
        (events.process_id == z_proc_id)
        & ((mu2tau_fakes_mask & (hcand_ele_mu_DM == -2)) | (e2tau_fakes_mask & (hcand_ele_mu_DM == -1)))
    )
    ztotautau = (
        (events.process_id == z_proc_id)
        & genuine_tau_mask
    )
    
    process_id = ak.where(ztoeeormm, ztoll_proc_id, events.process_id)
    process_id = ak.where(ztotautau, ztott_proc_id, events.process_id)
    
    #process_id = ak.where((mu2tau_fakes_mask & (hcand_ele_mu_DM == -2)), 51001, 51004)
    #process_id = ak.where((e2tau_fakes_mask & (hcand_ele_mu_DM == -1)), 51003, process_id)
    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int64)

    return events

