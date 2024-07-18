# coding: utf-8
"""
A producer to create Acoplanarity angles for differennt methods
"""

import os
from typing import Optional
from columnflow.production import Producer, producer
from columnflow.util import maybe_import

from httcp.production.PhiCP_Estimator import GetPhiCP
from columnflow.columnar_util import EMPTY_FLOAT, Route, set_ak_column, optional_column as optional

np = maybe_import("numpy")
ak = maybe_import("awkward")

#from IPython import embed



@producer(
    uses={
        "channel_id",
        "hcand.pt","hcand.eta", "hcand.phi", "hcand.mass", "hcand.decayMode", "hcand.charge",
        "hcand.IPx", "hcand.IPy", "hcand.IPz",
        "hcandprod.pt","hcandprod.eta", "hcandprod.phi", "hcandprod.mass", "hcandprod.pdgId",
        optional("GenTau.pt"), optional("GenTau.eta"), optional("GenTau.phi"), optional("GenTau.mass"),
        optional("GenTau.IPx"), optional("GenTau.IPy"), optional("GenTau.IPz"),
        optional("GenTau.decayMode"), optional("GenTau.charge"), 
        optional("GenTauProd.pt"), optional("GenTauProd.eta"), optional("GenTauProd.phi"), optional("GenTauProd.mass"),
    },
)
def ProducePhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> tuple[ak.Array,ak.Array,ak.Array]:

    PrepareP4 = lambda p4dict, mask : {key: ak.where(mask, val, val[:,:0]) for key, val in p4dict.items()}

    is_e      = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -1, axis=1), False)
    is_mu     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == -1, axis=1), False)
    is_pi     = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  0, axis=1), False)
    is_rho    = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  1, axis=1), False)
    is_a1DM2  = lambda leg: ak.fill_none(ak.firsts(leg.decayMode ==  2, axis=1), False)
    is_a1DM10 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 10, axis=1), False)
    is_a1DM11 = lambda leg: ak.fill_none(ak.firsts(leg.decayMode == 11, axis=1), False)
    is_a1     = is_a1DM10 # for now

    p4h1      = p4hcandinfo["p4h1"]
    p4h2      = p4hcandinfo["p4h2"]

    # etau
    mask_e_pi    = is_e(p4h1)   &   is_pi(p4h2)
    mask_e_rho   = is_e(p4h1)   &   is_rho(p4h2)
    mask_e_a1    = is_e(p4h1)   &   is_a1(p4h2)
    # mutau
    mask_mu_pi   = is_mu(p4h1)  &   is_pi(p4h2)
    mask_mu_rho  = is_mu(p4h1)  &   is_rho(p4h2)
    mask_mu_a1   = is_mu(p4h1)  &   is_a1(p4h2)
    # tautau
    mask_pi_pi   = is_pi(p4h1)  &   is_pi(p4h2)
    mask_pi_rho  = is_pi(p4h1)  &   is_rho(p4h2)
    mask_rho_pi  = is_rho(p4h1) &   is_pi(p4h2)
    mask_pi_a1   = is_pi(p4h1)  &   is_a1(p4h2)
    mask_a1_pi   = is_a1(p4h1)  &   is_pi(p4h2)
    mask_rho_rho = is_rho(p4h1) &   is_rho(p4h2)
    mask_rho_a1  = is_rho(p4h1) &   is_a1(p4h2)
    mask_a1_rho  = is_a1(p4h1)  &   is_rho(p4h2)
    mask_a1_a1   = is_a1(p4h1)  &   is_a1(p4h2)

    #print("BABUSHCHA")
    #embed()

    # GetPhiCP
    # Inputs:
    #   masked dict of p4hcand
    #   leg1_method, leg2_method
    #   leg1_mode,lmeg2_mode

    dummyPhiCP = ak.values_astype(events.event[:,None][:,:0], np.float32)

    # DPDP
    PhiCP_DPDP = ak.where(mask_rho_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_rho), "DP", "DP", "rho", "rho"), dummyPhiCP)
    PhiCP_DPDP = ak.where(mask_rho_a1,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1 ), "DP", "DP", "rho", "a1" ), PhiCP_DPDP)
    PhiCP_DPDP = ak.where(mask_a1_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1_rho ), "DP", "DP", "a1",  "rho"), PhiCP_DPDP)
    PhiCP_DPDP = ak.where(mask_a1_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_a1_a1  ), "DP", "DP", "a1",  "a1" ), PhiCP_DPDP)

    # PVPV
    PhiCP_PVPV = ak.where(mask_pi_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_pi  ), "PV", "PV", "pi",  "pi" ), dummyPhiCP)
    PhiCP_PVPV = ak.where(mask_pi_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho ), "PV", "PV", "pi",  "rho"), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_rho_pi,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_pi ), "PV", "PV", "rho", "pi" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_rho_rho, GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_rho), "PV", "PV", "rho", "rho"), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_pi_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1  ), "PV", "PV", "pi",  "a1" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_a1_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_a1_pi  ), "PV", "PV", "a1",  "pi" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_rho_a1,  GetPhiCP(PrepareP4(p4hcandinfo, mask_rho_a1 ), "PV", "PV", "rho", "a1" ), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_a1_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_a1_rho ), "PV", "PV", "a1",  "rho"), PhiCP_PVPV)
    PhiCP_PVPV = ak.where(mask_a1_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_a1_a1  ), "PV", "PV", "a1",  "a1" ), PhiCP_PVPV)

    # IPIP
    PhiCP_IPIP = ak.where(mask_e_pi,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_pi   ), "IP", "IP", "e",   "pi" ), dummyPhiCP)
    PhiCP_IPIP = ak.where(mask_mu_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_pi  ), "IP", "IP", "mu",  "pi" ), PhiCP_IPIP)
    PhiCP_IPIP = ak.where(mask_pi_pi,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_pi  ), "IP", "IP", "pi",  "pi" ), PhiCP_IPIP)

    # IPPV
    PhiCP_IPPV = ak.where(mask_e_rho,   GetPhiCP(PrepareP4(p4hcandinfo, mask_e_rho  ), "IP", "PV", "e",   "rho"), dummyPhiCP)
    PhiCP_IPPV = ak.where(mask_e_a1,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1   ), "IP", "PV", "e",   "a1" ), PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_mu_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_rho ), "IP", "PV", "mu",  "rho"), PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_mu_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1  ), "IP", "PV", "mu",  "a1" ), PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_pi_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho ), "IP", "PV", "pi",  "rho"), PhiCP_IPPV)
    PhiCP_IPPV = ak.where(mask_pi_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1  ), "IP", "PV", "pi",  "a1" ), PhiCP_IPPV)

    # IPDP
    PhiCP_IPDP = ak.where(mask_e_rho,   GetPhiCP(PrepareP4(p4hcandinfo, mask_e_rho  ), "IP", "DP", "e",   "rho"), dummyPhiCP)
    PhiCP_IPDP = ak.where(mask_e_a1,    GetPhiCP(PrepareP4(p4hcandinfo, mask_e_a1   ), "IP", "DP", "e",   "a1" ), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_mu_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_rho ), "IP", "DP", "mu",  "rho"), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_mu_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_mu_a1  ), "IP", "DP", "mu",  "a1" ), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_pi_rho,  GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_rho ), "IP", "DP", "pi",  "rho"), PhiCP_IPDP)
    PhiCP_IPDP = ak.where(mask_pi_a1,   GetPhiCP(PrepareP4(p4hcandinfo, mask_pi_a1  ), "IP", "DP", "pi",  "a1" ), PhiCP_IPDP)

    
    return events, PhiCP_PVPV, PhiCP_DPDP, PhiCP_IPIP, PhiCP_IPPV, PhiCP_IPDP



# ------------ DETECTOR LEVEL ----------- #
@producer(
    uses={
        ProducePhiCP,
    },
    produces={
        "PhiCP_PVPV", "PhiCP_DPDP", "PhiCP_IPIP", "PhiCP_IPPV", "PhiCP_IPDP",
    },
)
def ProduceDetPhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    events, PhiCP_PVPV, PhiCP_DPDP, PhiCP_IPIP, PhiCP_IPPV, PhiCP_IPDP = self[ProducePhiCP](events, p4hcandinfo)
    
    events = set_ak_column(events, "PhiCP_PVPV", PhiCP_PVPV)
    events = set_ak_column(events, "PhiCP_DPDP", PhiCP_DPDP)
    events = set_ak_column(events, "PhiCP_IPIP", PhiCP_IPIP)
    events = set_ak_column(events, "PhiCP_IPPV", PhiCP_IPPV)
    events = set_ak_column(events, "PhiCP_IPDP", PhiCP_IPDP)

    return events



# ------------ GENERATOR LEVEL ----------- #
@producer(
    uses={
        ProducePhiCP,
    },
    produces={
        "PhiCPGen_PVPV", "PhiCPGen_DPDP", "PhiCPGen_IPIP", "PhiCPGen_IPPV", "PhiCPGen_IPDP",
    },
    mc_only=True,
)
def ProduceGenPhiCP(
        self: Producer,
        events: ak.Array,
        p4hcandinfo: dict,
        **kwargs,
) -> ak.Array:
    events, PhiCP_PVPV, PhiCP_DPDP, PhiCP_IPIP, PhiCP_IPPV, PhiCP_IPDP = self[ProducePhiCP](events, p4hcandinfo)
    PhiCP_IPIP = ak.nan_to_num(PhiCP_IPIP, nan=0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPPV = ak.nan_to_num(PhiCP_IPPV, nan=0.0) # WRONG: CHECK IP FOR GEN PARTICLES
    PhiCP_IPDP = ak.nan_to_num(PhiCP_IPDP, nan=0.0) # WRONG: CHECK IP FOR GEN PARTICLES

    events = set_ak_column(events, "PhiCPGen_PVPV", PhiCP_PVPV)
    events = set_ak_column(events, "PhiCPGen_DPDP", PhiCP_DPDP)
    events = set_ak_column(events, "PhiCPGen_IPIP", PhiCP_IPIP)
    events = set_ak_column(events, "PhiCPGen_IPPV", PhiCP_IPPV)
    events = set_ak_column(events, "PhiCPGen_IPDP", PhiCP_IPDP)

    return events
