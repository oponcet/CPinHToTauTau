# coding: utf-8

"""
main calibration script
"""

from columnflow.calibration import Calibrator, calibrator
from columnflow.calibration.cms.jets import jec, jec_nominal, jer
from columnflow.production.cms.seeds import deterministic_seeds
from columnflow.util import maybe_import
from columnflow.columnar_util import set_ak_column

from httcp.calibration.tau import tau_energy_scale
from httcp.util import IF_RUN2, IF_RUN3

np = maybe_import("numpy")
ak = maybe_import("awkward")

# derive calibrators to add settings
#jec_full = jec.derive("jec_full", cls_dict={"mc_only": True, "nominal_only": True})

@calibrator(
    uses={
        deterministic_seeds,
        #jec_nominal, jec_full, jer,
        tau_energy_scale,
    },
    produces={
        deterministic_seeds,
        #jec_nominal, jec_full, jer,
        tau_energy_scale,
    },
)
def main(self: Calibrator, events: ak.Array, **kwargs) -> ak.Array:
    events = self[deterministic_seeds](events, **kwargs)

    #if self.dataset_inst.is_data or not self.global_shift_inst.is_nominal:
    #    events = self[jec_nominal](events, **kwargs)
    #else:
    #    events = self[jec_full](events, **kwargs)
    #    events = self[jer](events, **kwargs)

    if self.dataset_inst.is_mc: 
        ##Apply tau energy scale correction
        events = self[tau_energy_scale](events, **kwargs)


    return events
