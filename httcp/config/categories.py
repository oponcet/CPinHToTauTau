
# coding: utf-8

# To handle many categories:
# https://github.com/columnflow/columnflow/commit/3a104d633fa47a8efc789f7aba054ed967017347#diff-01da7ecfbc4b8bb83460821147201da604f9825e32cbd51e365bdfcbc7cb0912R61
# https://github.com/columnflow/columnflow/issues/547


import law
import order as od

from columnflow.config_util import add_category, create_category_combinations
from httcp.util import call_once_on_config

logger = law.logger.get_logger(__name__)


@call_once_on_config()
def add_leptauh_ABCD_categories(config: od.Config) -> None:
    add_category(config,name="lepA",  id=101, selection="cat_ss_iso2_bveto_lowmt",      label="A")
    add_category(config,name="lepB",  id=102, selection="cat_ss_noniso2_bveto_lowmt",   label="B")
    add_category(config,name="lepA0", id=103, selection="cat_os_iso2_nobveto_lowmt",    label="A0")
    add_category(config,name="lepB0", id=104, selection="cat_os_noniso2_nobveto_lowmt", label="B0")
    add_category(config,name="lepA1", id=105, selection="cat_os_iso2_bveto_highmt",     label="A1")
    add_category(config,name="lepB1", id=106, selection="cat_os_iso2_bveto_lowmt",      label="B1")
    add_category(config,name="lepD",  id=107, selection="cat_os_iso2_bveto_lowmt",      label="D")
    add_category(config,name="lepC",  id=108, selection="cat_os_noniso2_bveto_lowmt",   label="C")

@call_once_on_config()
def add_lept_cp_categories(config: od.Config) -> None:
    add_category(config, name="pi",        id=201, selection="cat_pi",     label=r"$\pi$")  # h2 -> pi
    add_category(config, name="rho",       id=202, selection="cat_rho",    label=r"$\rho$") # h2 -> rho
    add_category(config, name="a1dm2",     id=203, selection="cat_a1dm2",  label=r"$a_{1}(2p-2\pi^{0})$") # h2 -> a1
    add_category(config, name="a1dm10",    id=204, selection="cat_a1dm10", label=r"$a_{1}(3p-0\pi^{0})$") # h2 -> a1
    add_category(config, name="a1dm11",    id=205, selection="cat_a1dm11", label=r"$a_{1}(3p-1\pi^{0})$") # h2 -> a1

@call_once_on_config()
def add_hadtauh_ABCD_categories(config: od.Config) -> None:
    add_category(config,name="hadA",  id=301, selection="cat_ss_iso1_iso2_bveto",      label="A")
    add_category(config,name="hadB",  id=302, selection="cat_ss_noniso1_iso2_bveto",   label="B")
    add_category(config,name="hadA0", id=303, selection="cat_ss_iso1_noniso2_bveto",   label="A0")
    add_category(config,name="hadB0", id=304, selection="cat_ss_noniso1_noniso2_bveto",label="B0")
    add_category(config,name="hadD0", id=305, selection="cat_os_iso1_noniso2_bveto",   label="D0")
    add_category(config,name="hadC0", id=306, selection="cat_os_noniso1_noniso2_bveto",label="C0")
    add_category(config,name="hadD",  id=307, selection="cat_os_iso1_iso2_bveto",      label="D")
    add_category(config,name="hadC",  id=308, selection="cat_os_noniso1_iso2_bveto",   label="C")
    
@call_once_on_config()
def add_hadr_cp_categories(config: od.Config) -> None:
    add_category(config, name="pi_pi",         id=401, selection="cat_pi_pi",           label=r"$\pi-\pi$")  
    add_category(config, name="pi_rho",        id=402, selection="cat_pi_rho",          label=r"$\pi-\rho$") # h1 -> rho
    add_category(config, name="pi_a1dm2",      id=403, selection="cat_pi_a1dm2",        label=r"$\pi-a_{1}(2p2\pi^{0})$") # h1 -> a1
    add_category(config, name="pi_a1dm10",     id=404, selection="cat_pi_a1dm10",       label=r"$\pi-a_{1}(3p0\pi^{0})$") # h1 -> a1
    add_category(config, name="pi_a1dm11",     id=405, selection="cat_pi_a1dm11",       label=r"$\pi-a_{1}(3p1\pi^{0})$") # h1 -> a1
    add_category(config, name="rho_rho",       id=406, selection="cat_rho_rho",         label=r"$\rho-\rho$") # h1 -> rho
    add_category(config, name="rho_a1dm2",     id=407, selection="cat_rho_a1dm2",       label=r"$\rho-a_{1}(2p2\pi^{0})$") # h1 -> a1
    add_category(config, name="rho_a1dm10",    id=408, selection="cat_rho_a1dm10",      label=r"$\rho-a_{1}(3p0\pi^{0})$") # h1 -> a1
    add_category(config, name="rho_a1dm11",    id=409, selection="cat_rho_a1dm11",      label=r"$\rho-a_{1}(3p1\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm2_a1dm2",   id=410, selection="cat_a1dm2_a1dm2",     label=r"$a_{1}-a_{1}(2p2\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm2_a1dm10",  id=411, selection="cat_a1dm2_a1dm10",    label=r"$a_{1}(2p2\pi^{0})-a_{1}(3p0\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm2_a1dm11",  id=412, selection="cat_a1dm2_a1dm11",    label=r"$a_{1}(2p2\pi^{0})-a_{1}(3p1\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm10_a1dm10", id=413, selection="cat_a1dm10_a1dm10",   label=r"$a_{1}-a_{1}(3p0\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm10_a1dm11", id=414, selection="cat_a1dm10_a1dm11",   label=r"$a_{1}(3p0\pi^{0})-a_{1}(3p1\pi^{0})$") # h1 -> a1
    add_category(config, name="a1dm11_a1dm11", id=415, selection="cat_a1dm11_a1dm11",   label=r"$a_{1}(3p1\pi^{0})-a_{1}(3p1\pi^{0})$") # h1 -> a1
    
    


@call_once_on_config()
def add_etau_mutau_categories(config: od.Config) -> None:
    add_category(config,
                 name="etau",
                 id=1000,
                 selection="cat_etau",
                 label=r"$e\tau_{h}$")
    add_category(config,
                 name="mutau",
                 id=2000,
                 selection="cat_mutau",
                 label=r"$\mu\tau_{h}$")

    add_leptauh_ABCD_categories(config)
    add_lept_cp_categories(config)

    categories = {
        "channel": [config.get_category("etau"),
                    config.get_category("mutau")],
        "abcd"   : [config.get_category("lepA"),  config.get_category("lepB"),
                    config.get_category("lepA0"), config.get_category("lepB0"),
                    config.get_category("lepA1"), config.get_category("lepB1"),
                    config.get_category("lepC"),  config.get_category("lepD")],
        #"abcd"   : [config.get_category("lepC"),  config.get_category("lepD")],
        #"cp"     : [config.get_category("pi"),
        #            config.get_category("rho"),
        #            config.get_category("a1dm2"),
        #            config.get_category("a1dm10"),
        #            config.get_category("a1dm11")],
    }
    
    def name_fn(categories):
        catlist = [cat.name for cat in categories.values() if cat]
        catname = "__".join(catlist)
        return catname

    
    def kwargs_fn(categories):
        dictnry = {
            # just increment the category id
            # NOTE: for this to be deterministic, the order of the categories must no change!
            #"id": "+",
            "id": sum([c.id for c in categories.values()]),
            # join all tags
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
            # auxiliary information
            #"aux": {
            #    # the fake factors group name
            #    "fakefactors_group": name_fn({name: cat for name, cat in categories.items() if name not in {"sign", "tau2"}}),
            #},
        }
        return dictnry
    
    create_category_combinations(config, categories, name_fn, kwargs_fn)
    
    
    
@call_once_on_config()
def add_tautau_categories(config: od.Config) -> None:
    add_category(config,
                 name="tautau",
                 id=4000,
                 selection="cat_tautau",
                 label=r"$\tau_{h}\tau_{h}$")

    add_hadtauh_ABCD_categories(config)
    add_hadr_cp_categories(config)
    
    categories = {
        "channel": [config.get_category("tautau")],
        "abcd"   : [config.get_category("hadA"),  config.get_category("hadB"),
                    config.get_category("hadA0"), config.get_category("hadB0"),
                    config.get_category("hadC0"), config.get_category("hadD0"),
                    config.get_category("hadC"),  config.get_category("hadD")],
        #"abcd"   : [config.get_category("hadC"),  config.get_category("hadD")],
        #"cp"     : [config.get_category("pi_pi"),
        #            config.get_category("pi_rho"),
        #            config.get_category("pi_a1dm2"),
        #            config.get_category("pi_a1dm10"),
        #            config.get_category("pi_a1dm11"),
        #            config.get_category("rho_rho"),
        #            config.get_category("rho_a1dm2"),
        #            config.get_category("rho_a1dm10"),
        #            config.get_category("rho_a1dm11"),
        #            config.get_category("a1dm2_a1dm10"),
        #            config.get_category("a1dm2_a1dm11"),
        #            config.get_category("a1dm10_a1dm10"),
        #            config.get_category("a1dm10_a1dm11"),
        #            config.get_category("a1dm11_a1dm11")],
    }
    
    def name_fn(categories):
        catlist = [cat.name for cat in categories.values() if cat]
        catname = "__".join(catlist)
        #print(catname)
        return catname
    
    def kwargs_fn(categories):
        dictnry = {
            # just increment the category id
            # NOTE: for this to be deterministic, the order of the categories must no change!
            #"id": "+",
            "id": sum([c.id for c in categories.values()]),
            # join all tags
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
            # auxiliary information
            #"aux": {
            #    # the fake factors group name
            #    "fakefactors_group": name_fn({name: cat for name, cat in categories.items() if name not in {"sign", "tau2"}}),
            #},
        }
        #print("ID : ", dictnry["id"])
        return dictnry
        
    create_category_combinations(config, categories, name_fn, kwargs_fn)
    
    
@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_category(config, name="incl", id=9, selection="cat_incl", label=r"$Inclusive$")

    add_etau_mutau_categories(config)
    add_tautau_categories(config)
