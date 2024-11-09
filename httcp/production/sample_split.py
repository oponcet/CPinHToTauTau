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
        "unknown"   : 0,
        "prompt_e"  : 1,
        "prompt_mu" : 2,
        "tau->e"    : 3,
        "tau->mu"   : 4,
        "tau_had"   : 5,
        "fake"      : 6
    }
    #hcand_ele_mu_DM = events.hcand.decayMode[:,0]
    hcand_tau_DM = events.hcand.decayMode[:,1]
    
    # hcand_Tau_idx = events.hcand.rawIdx[:,1:2]
    # mask_hcand_tau_idx = ak.flatten(hcand_Tau_idx)
    match = events.Tau.genPartFlav
    #### The ak.any() is a nasty fix that need to be removed, we have tautau that needs to be taken into account properly
    mu2tau_fakes_mask = ak.any(((match == tau_part_flav["prompt_mu"]) | (match == tau_part_flav["tau->mu"])),axis=1)
    e2tau_fakes_mask  = ak.any(((match == tau_part_flav["prompt_e"]) | (match == tau_part_flav["tau->e"])),axis=1)
    j2tau_fakes_mask  = ak.any(((match == tau_part_flav["unknown"]) | (match == tau_part_flav["fake"])),axis=1)
    genuine_tau_mask  = ak.any((match == tau_part_flav["tau_had"]),axis=1)
    
    z_proc_id       = 51100 #self.dataset_inst.processes.values()[0].id

    # N.B. For tau-tau channel, lets assume the leading tau is true
    z_j2tau_proc_id = 51097 # if the hcand2 i.e. tauh is from jet or not
    z_t2tau_proc_id = 51098 # if the hcand2 i.e. tauh is true or not
    z_l2tau_proc_id = 51099 # if the hcand2 i.e. tauh is from e/mu or not

    #ztoee_proc_id = 51101
    #ztomm_proc_id = 51102
    #ztounknown_id = 51097
    #ztott_proc_id = 51098
    #ztoll_proc_id = 51099
    
    #from IPython import embed; embed()
    
    #process_id = np.array(events.process_id, dtype=np.int64)

    #ztoeeormm = (
    #    (events.process_id == z_proc_id)
    #    & ((mu2tau_fakes_mask & (hcand_ele_mu_DM == -2)) | (e2tau_fakes_mask & (hcand_ele_mu_DM == -1)))
    #)
    #ztotautau = (
    #    (events.process_id == z_proc_id)
    #    & genuine_tau_mask
    #)
    
    #process_id = ak.where((mu2tau_fakes_mask | e2tau_fakes_mask), ztoll_proc_id, events.process_id)
    #process_id = ak.where(genuine_tau_mask, ztott_proc_id, process_id)
    #process_id = ak.where(unknown_mask, ztounknown_id, process_id)

    process_id = ak.where((mu2tau_fakes_mask | e2tau_fakes_mask),
                          z_l2tau_proc_id,
                          events.process_id)
    process_id = ak.where(j2tau_fakes_mask,
                          z_j2tau_proc_id,
                          process_id)
    process_id = ak.where(genuine_tau_mask,
                          z_t2tau_proc_id,
                          process_id)

    events = remove_ak_column(events, "process_id")
    events = set_ak_column(events, "process_id", process_id, value_type=np.int64)

    return events
