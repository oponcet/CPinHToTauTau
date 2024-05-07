import os
import sys
from typing import Optional
from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")
maybe_import("coffea.nanoevents.methods.nanoaod")

def convert_to_coffea_p4(zipped_item):
    return ak.zip(
        zipped_item,
        with_name = "PtEtaPhiMLorentzVector",
        behavior  = coffea.nanoevents.methods.vector.behavior,
    )

def reconstructPi0(photons, method: Optional[str] = "simpleIC"):
    photons = ak.with_name(photons, "PtEtaPhiMLorentzVector")
    photons_sorted_pt_indices = ak.argsort(photons.pt, ascending=False)
    photons = photons[photons_sorted_pt_indices]

    p4_pi0 = None

    if method == "simpleIC":
        photons_px = ak.sum(photons.px, axis=1)
        photons_py = ak.sum(photons.py, axis=1)
        photons_pt = np.sqrt(photons_px ** 2 + photons_py ** 2)

        #pt_pi0    = photons[:, 0:1].pt
        pt_pi0    = photons_pt
        eta_pi0   = photons[:, 0:1].eta
        phi_pi0   = photons[:, 0:1].phi
        pdgid_pi0 = ak.values_astype(111 * ak.ones_like(eta_pi0), "int64")
        mass_pi0  = 0.135 * ak.ones_like(eta_pi0)
        
        p4_pi0 = convert_to_coffea_p4({
            "pt"    : pt_pi0,
            "eta"   : eta_pi0,
            "phi"   : phi_pi0,
            "mass"  : mass_pi0,
            "pdgId" : pdgid_pi0,
        })
        
    elif method == "simpleMB":
        has_one_photon = ak.num(photons.pt, axis=1) == 1
        pass

    return p4_pi0
