# coding: utf-8

"""
Exemplary selection methods.
"""

from columnflow.categorization import Categorizer, categorizer
from columnflow.util import maybe_import

ak = maybe_import("awkward")


#
# categorizer functions used by categories definitions
#

@categorizer(uses={"event"})
def cat_incl(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # fully inclusive selection
    return events, ak.ones_like(events.event) == 1


@categorizer(uses={"Jet.pt"})
def cat_2j(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # two or more jets
    return events, ak.num(events.Jet.pt, axis=1) >= 2

@categorizer(uses={"channel_id"})
def sel_etau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # ee channel
    ch = self.config_inst.get_channel("etau")
    return events, events["channel_id"] == ch.id


@categorizer(uses={"channel_id"})
def sel_mutau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # mm channel
    ch = self.config_inst.get_channel("mutau")
    return events, events["channel_id"] == ch.id

@categorizer(uses={"channel_id"})
def sel_tautau(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    # mm channel
    ch = self.config_inst.get_channel("tautau")
    return events, events["channel_id"] == ch.id

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_c(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Control region ( iso < 0.15, same sign pair)
    sel = (events.rel_charge > 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) < 0.15)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_d(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Signal region ( iso < 0.15, opposite sign pair)
    sel = (events.rel_charge < 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) < 0.15)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_a(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Region for transfer factor calculation( iso > 0.15, same sign pair)
    sel = (events.rel_charge > 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) >= 0.15)  \
        & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) <= 0.30)
    return events, sel

@categorizer(uses={"rel_charge", "Muon.pfRelIso04_all"})
def cat_b(self: Categorizer, events: ak.Array, **kwargs) -> tuple[ak.Array, ak.Array]:
    #Region for transfer factor calculation( iso > 0.15, opposite sign pair)
    sel = (events.rel_charge < 0) & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) >= 0.15) \
        & (ak.firsts(events.Muon.pfRelIso04_all, axis=1) <= 0.30)
    return events, sel
