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

def getMaxEtaTauStrip(pt):
    temp = 0.20 * np.power(pt, -0.66)
    ref1 = ak.where(temp > 0.15, 0.15, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2

def getMaxPhiTauStrip(pt):
    temp = 0.35 * np.power(pt, -0.71)
    ref1 = ak.where(temp > 0.30, 0.30, temp)
    ref2 = ak.where(ref1 > 0.05, ref1, 0.05)
    return ref2


def reconstructPi0(hcandp4, photons, method: Optional[str] = "simpleIC"):
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
        has_atleast_one_photon = ak.num(photons.pt, axis=1) > 0
        hcandp4 = ak.where(has_atleast_one_photon, hcandp4, hcandp4[:,:0])
        photons_p4 = ak.where(has_atleast_one_photon,
                              photons_p4[:,0:1],
                              photons_p4[:,:0])
        
        deta_photons_hcand = (photons_p4).metric_table(hcandp4, metric = lambda a,b: np.abs(a.eta - b.eta))
        dphi_photons_hcand = (photons_p4).metric_table(hcandp4, metric = lambda a,b: np.abs(a.delta_phi(b)))
        
        maxeta_photons = getMaxEtaTauStrip(photons_p4.pt)
        maxphi_photons = getMaxPhiTauStrip(photons_p4.pt)

        mask_photons = ((np.abs(deta_photons_hcand) < maxeta_photons)
                        & (np.abs(dphi_photons_hcand) < maxphi_photons))
        
        strip_photons_p4 = convert_to_coffea_p4(
            {
                "pt"   : photons_p4.pt[mask_photons],
                "eta"  : photons_p4.eta[mask_photons],
                "phi"  : photons_p4.phi[mask_photons],
                "mass" : photons_p4.mass[mask_photons],
            }
        )

        has_one_photon = ak.num(strip_photons_p4.pt, axis=1) == 1
        
        strip_photons_p4_pair = ak.combinations(strip_photons_p4, 2, axis=1)
        strip_photons_p4_pair_0, strip_photons_p4_pair_1 = ak.unzip(strip_photons_p4_pair)
        strip_photons_mass = (strip_photons_p4_pair_0 + strip_photons_p4_pair_1).mass
        mass_val = np.abs(strip_photons_mass - pi0RecoM)
        mass_sorted_idx = ak.argsort(mass_val, axis=1)
        strip_photons_p4_pair_sorted = strip_photons_p4_pair[mass_sorted_idx]
        mass_mask = mass_val < 2 * pi0RecoW
        evt_mask_no_pair = ak.sum(mass_mask, axis=1) == 0
        strip_photons_p4_pair_sorted_pass_mass = strip_photons_p4_pair_sorted[mass_mask]
        strip_photons_p4_mass_selected = ak.concatenate([strip_photons_p4_pair_sorted_pass_mass["0"][:,:1],
                                                         strip_photons_p4_pair_sorted_pass_mass["1"][:,:1]], axis=1)

        sel_strip_photons_p4 = ak.where(has_one_photon,
                                        strip_photons_p4,
                                        ak.where(evt_mask_no_pair,
                                                 strip_photons_p4[:,:1],
                                                 strip_photons_p4_mass_selected)
                                    )
        
        sel_strip_pizero_p4 = ak.from_regular(ak.sum(sel_strip_photons_p4, axis=-1)[:,None])

        p4_pi0 = sel_strip_pizero_p4


    return p4_pi0
