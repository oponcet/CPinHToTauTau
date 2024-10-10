
# coding: utf-8

"""
Definition of categories.
- etau__os__iso
- etau__os__noniso
- etau__incl__ss__iso
- etau__incl__ss__noniso
- etau__2j__os__iso
- etau__2j__os__noniso
- etau__2j__ss__iso
- etau__2j__ss__noniso
- mutau__incl__os__iso
- mutau__incl__os__noniso
- mutau__incl__ss__iso
- mutau__incl__ss__noniso
- mutau__2j__os__iso
- mutau__2j__os__noniso
- mutau__2j__ss__iso
- mutau__2j__ss__noniso
- tautau__os__iso : tautau incl opposite sign both tau passing MediumvsJet WP = SR (D)
- tautau__incl__os__noniso : tautau incl opposite sign both tau passing VVLoosevsJet WP but not Medium = AR (C)
- tautau__incl__ss__iso : tautau incl same sign both tau passing MediumvsJet WP = DR (A)
- tautau__incl__ss__noniso : tautau incl same sign both tau passing VVLoosevsJet WP but not Medium = DR (B)
- tautau__2j__os__iso
- tautau__2j__os__noniso
- tautau__2j__ss__iso
- tautau__2j__ss__noniso

"""

import order as od

from columnflow.config_util import add_category, create_category_combinations


def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    # lepton channels
    add_category(config, name="etau", id=1, selection="sel_etau", label=r"$e\tau_{h}$")
    add_category(config, name="mutau", id=2, selection="sel_mutau", label=r"$\mu\tau_{h}$")
    add_category(config, name="tautau", id=4, selection="sel_tautau", label=r"$\tau_{h}\tau_{h}$")

    # Fake factors regions
    add_category(config, name="os",      id=10, selection="cat_os",       label="Opposite sign", tags={"os"})
    add_category(config, name="ss",      id=11, selection="cat_ss",       label="Same sign", tags={"ss"})
    add_category(config, name="iso1",    id=12, selection="cat_iso_1",    label=r"$\tau_{h,1} isolated$", tags={"iso1"})
    add_category(config, name="noniso1", id=13, selection="cat_noniso_1", label=r"$\tau_{h,1} non-isolated$", tags={"noniso1"})  # noqa: E501
    add_category(config, name="iso2",    id=14, selection="cat_iso_2",    label=r"$\tau_{h,2} isolated$", tags={"iso2"})
    add_category(config, name="noniso2", id=15, selection="cat_noniso_2", label=r"$\tau_{h,2} non-isolated$", tags={"noniso2"})  # noqa: E501
    add_category(config, name="lowmt",   id=16, selection="cat_low_mt",   label=r"$low m_{T} (< 50 GeV)$", tags={"lowmt"})
    add_category(config, name="highmt",  id=17, selection="cat_high_mt",  label=r"$high m_{T} (> 50 GeV)$", tags={"highmt"})
    add_category(config, name="nobjet",  id=18, selection="cat_has_no_b", label=r"$0 b-jet$", tags={"nobjet"})
    add_category(config, name="bjet",    id=19, selection="cat_has_b",    label=r"$> 0 b-jet$", tags={"bjet"})

    
    # kinematic categories
    #add_category(config, name="incl", id=100, selection="cat_incl", label="inclusive")
    #add_category(config, name="2j", id=110, selection="cat_2j", label="2 jets")

    #
    # build groups
    #

    categories = {
        # channels first
        "channel": [config.get_category("etau"), config.get_category("mutau"), config.get_category("tautau")],
        # kinematic regions in the middle (to be extended)
        #"kin": [config.get_category("incl"), config.get_category("2j")],
        # Fake factors regions last
        "sign": [config.get_category("os"), config.get_category("ss")],
        "tau1": [config.get_category("iso1"), config.get_category("noniso1")],
        "tau2": [config.get_category("iso2"), config.get_category("noniso2")],
        "b"   : [config.get_category("bjet"), config.get_category("nobjet")],
    }

    def name_fn(categories):
        return "__".join(cat.name for cat in categories.values() if cat)

    def kwargs_fn(categories):
        return {
            # just increment the category id
            # NOTE: for this to be deterministic, the order of the categories must no change!
            "id": "+",
            # join all tags
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
            # auxiliary information
            #"aux": {
            #    # the fake factors group name
            #    "fakefactors_group": name_fn({name: cat for name, cat in categories.items() if name not in {"sign", "tau2"}}),
            #},
        }

    create_category_combinations(config, categories, name_fn, kwargs_fn)
