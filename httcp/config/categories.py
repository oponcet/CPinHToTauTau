# coding: utf-8

# To handle many categories:
# https://github.com/columnflow/columnflow/commit/3a104d633fa47a8efc789f7aba054ed967017347#diff-01da7ecfbc4b8bb83460821147201da604f9825e32cbd51e365bdfcbc7cb0912R61
# https://github.com/columnflow/columnflow/issues/547
# https://github.com/columnflow/columnflow/issues/559

import law
import order as od

from columnflow.config_util import add_category, create_category_combinations
from httcp.util import call_once_on_config

logger = law.logger.get_logger(__name__)


# ############################################################### #
# To create combinations of categories                            #
# the entire root -> leaf categories will be in the category_ids  #
# other combinatorics will be produced in createHistograms task,  #
# if mentioned                                                    #
# ############################################################### #
def create_nested_categories(config: od.Config, categories) -> None:

    def name_fn(categories):
        catlist = [cat.name for cat in categories.values() if cat]
        catname = "__".join(catlist)
        return catname
    
    def kwargs_fn(categories):
        return {
            "id": sum([c.id for c in categories.values()]),
            "tags": set.union(*[cat.tags for cat in categories.values() if cat]),
        }

    n = create_category_combinations(config, categories, name_fn, kwargs_fn)
    logger.info(f"\t{n} categories created")



# ###################### #
#   for tau-tau channel  #
# ###################### #
@call_once_on_config()
def add_hadtauh_ABCD_categories(config: od.Config) -> None:
    add_category(config,name="hadA",  id=100,  selection="cat_ss_iso1_iso2_bveto",       label="A",  tags={"ss","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadB",  id=300,  selection="cat_ss_noniso1_iso2_bveto",    label="B",  tags={"ss","noniso1","iso2",   "bveto"})
    add_category(config,name="hadA0", id=500,  selection="cat_ss_iso1_noniso2_bveto",    label="A0", tags={"ss","iso1",   "noniso2","bveto"})
    add_category(config,name="hadB0", id=700,  selection="cat_ss_noniso1_noniso2_bveto", label="B0", tags={"ss","noniso1","noniso2","bveto"})
    add_category(config,name="hadD0", id=900,  selection="cat_os_iso1_noniso2_bveto",    label="D0", tags={"os","iso1",   "noniso2","bveto"})
    add_category(config,name="hadC0", id=1100, selection="cat_os_noniso1_noniso2_bveto", label="C0", tags={"os","noniso1","noniso2","bveto"})
    add_category(config,name="hadD",  id=1300, selection="cat_os_iso1_iso2_bveto",       label="D",  tags={"os","iso1",   "iso2",   "bveto"})
    add_category(config,name="hadC",  id=1500, selection="cat_os_noniso1_iso2_bveto",    label="C",  tags={"os","noniso1","iso2",   "bveto"})
    
@call_once_on_config()
def add_hadr_cp_categories(config: od.Config) -> None:
    # add_category(config, name="pi_1",          id=1,  selection="cat_pi_1",            label=r"$\pi$",                                 tags={"tau1pi"         })  # h1 -> pi
    # add_category(config, name="rho_1",         id=3,  selection="cat_rho_1",           label=r"$\rho$",                                tags={"tau1rho"        })  # h1 -> rho
    # add_category(config, name="a1dm2_1",       id=5,  selection="cat_a1dm2_1",         label=r"$a_{1}(2p-2\pi^{0})$",                  tags={"tau1a1DM2"      })  # h1 -> a1
    # add_category(config, name="a1dm10_1",      id=7,  selection="cat_a1dm10_1",        label=r"$a_{1}(3p-0\pi^{0})$",                  tags={"tau1a1DM10"     })  # h1 -> a1
    # add_category(config, name="a1dm11_1",      id=9,  selection="cat_a1dm11_1",        label=r"$a_{1}(3p-1\pi^{0})$",                  tags={"tau1a1DM11"     })  # h1 -> a1
    # add_category(config, name="pi_rho_1",      id=11, selection="cat_pi_rho_1",        label=r"$\pi/\rho$",                            tags={"tau1piorrho"    })  # h1 -> a1
    add_category(config, name="pi_pi",         id=11, selection="cat_pi_pi",           label=r"$\pi-\pi$",                             tags={"pi","pi"        })  # 
    add_category(config, name="pi_rho",        id=13, selection="cat_pi_rho",          label=r"$\pi-\rho$",                            tags={"pi","rho"       })  # 
    # add_category(config, name="pi_a1dm2",      id=15, selection="cat_pi_a1dm2",        label=r"$\pi-a_{1}(2p2\pi^{0})$",               tags={"pi","a1DM2"     })  # 
    add_category(config, name="pi_a1dm10",     id=17, selection="cat_pi_a1dm10",       label=r"$\pi-a_{1}(3p0\pi^{0})$",               tags={"pi","a1DM10"    })  # 
    add_category(config, name="pi_a1dm11",     id=19, selection="cat_pi_a1dm11",       label=r"$\pi-a_{1}(3p1\pi^{0})$",               tags={"pi","a1DM11"    })  # 
    add_category(config, name="rho_rho",       id=21, selection="cat_rho_rho",         label=r"$\rho-\rho$",                           tags={"rho","rho"      })  # 
    # add_category(config, name="rho_a1dm2",     id=23, selection="cat_rho_a1dm2",       label=r"$\rho-a_{1}(2p2\pi^{0})$",              tags={"rho","a1DM2"    })  # 
    add_category(config, name="rho_a1dm10",    id=25, selection="cat_rho_a1dm10",      label=r"$\rho-a_{1}(3p0\pi^{0})$",              tags={"rho","a1DM10"   })  # 
    add_category(config, name="rho_a1dm11",    id=27, selection="cat_rho_a1dm11",      label=r"$\rho-a_{1}(3p1\pi^{0})$",              tags={"rho","a1DM11"   })  # 
    # add_category(config, name="a1dm2_a1dm2",   id=29, selection="cat_a1dm2_a1dm2",     label=r"$a_{1}-a_{1}(2p2\pi^{0})$",             tags={"a1DM2","a1DM2"  })  # 
    # add_category(config, name="a1dm2_a1dm10",  id=31, selection="cat_a1dm2_a1dm10",    label=r"$a_{1}(2p2\pi^{0})-a_{1}(3p0\pi^{0})$", tags={"a1DM2","a1DM10" })  # 
    # add_category(config, name="a1dm2_a1dm11",  id=33, selection="cat_a1dm2_a1dm11",    label=r"$a_{1}(2p2\pi^{0})-a_{1}(3p1\pi^{0})$", tags={"a1DM2","a1DM11" })  # 
    add_category(config, name="a1dm10_a1dm10", id=35, selection="cat_a1dm10_a1dm10",   label=r"$a_{1}-a_{1}(3p0\pi^{0})$",             tags={"a1DM10","a1DM10"})  # 
    add_category(config, name="a1dm10_a1dm11", id=37, selection="cat_a1dm10_a1dm11",   label=r"$a_{1}(3p0\pi^{0})-a_{1}(3p1\pi^{0})$", tags={"a1DM10","a1DM11"})  # 
    add_category(config, name="a1dm11_a1dm11", id=39, selection="cat_a1dm11_a1dm11",   label=r"$a_{1}(3p1\pi^{0})-a_{1}(3p1\pi^{0})$", tags={"a1DM11","a1DM11"})  # 


@call_once_on_config()
def add_tautau_real_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("tautau")],
        "TorF"   : [config.get_category("real_1")],
        "abcd"   : [
            config.get_category("hadA"),  config.get_category("hadB"),
            config.get_category("hadA0"), config.get_category("hadB0"),
            config.get_category("hadC0"), config.get_category("hadD0"),
            config.get_category("hadC"),  config.get_category("hadD"),
        ],
        # "nj"     : [
        #     config.get_category("has_0j"), config.get_category("has_1j"), config.get_category("has_2j"),
        # ],
        "cp"     : [
            # config.get_category("pi_1"),
            # config.get_category("rho_1"),
            # config.get_category("a1dm2_1"),
            # config.get_category("a1dm10_1"),
            # config.get_category("a1dm11_1"),
            config.get_category("pi_pi"),
            config.get_category("pi_rho"),
            # config.get_category("pi_a1dm2"),
            config.get_category("pi_a1dm10"),
            config.get_category("pi_a1dm11"),
            config.get_category("rho_rho"),
            # config.get_category("rho_a1dm2"),
            config.get_category("rho_a1dm10"),
            config.get_category("rho_a1dm11"),
            # config.get_category("a1dm2_a1dm10"),
            # config.get_category("a1dm2_a1dm11"),
            config.get_category("a1dm10_a1dm10"),
            config.get_category("a1dm10_a1dm11"),
            config.get_category("a1dm11_a1dm11"),
        ],
    }
    logger.info("tautau_real_categories")
    create_nested_categories(config, categories)
    

@call_once_on_config()
def add_tautau_fake_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("tautau")],
        "TorF"   : [config.get_category("fake_1")],
        "abcd"   : [
            config.get_category("hadA"),  config.get_category("hadB"),
            config.get_category("hadA0"), config.get_category("hadB0"),
            config.get_category("hadC0"), config.get_category("hadD0"),
            config.get_category("hadC"),  config.get_category("hadD"),
        ],
        "nj"     : [
            config.get_category("has_0j"), config.get_category("has_1j"), config.get_category("has_2j"),
        ],
        "cp"     : [
            config.get_category("pi_1"),
            config.get_category("rho_1"),
            config.get_category("a1dm2_1"),
            config.get_category("a1dm10_1"),
            config.get_category("a1dm11_1"),
            config.get_category("pi_pi"),
            config.get_category("pi_rho"),
            config.get_category("pi_a1dm2"),
            config.get_category("pi_a1dm10"),
            config.get_category("pi_a1dm11"),
            #config.get_category("rho_rho"),
            #config.get_category("rho_a1dm2"),
            #config.get_category("rho_a1dm10"),
            #config.get_category("rho_a1dm11"),
            #config.get_category("a1dm2_a1dm10"),
            #config.get_category("a1dm2_a1dm11"),
            #config.get_category("a1dm10_a1dm10"),
            #config.get_category("a1dm10_a1dm11"),
            #config.get_category("a1dm11_a1dm11"),
        ],
    }
    logger.info("tautau_fake_categories")
    create_nested_categories(config, categories)



    
# ###################### #
#   for lep-tau channel  #
# ###################### #
@call_once_on_config()
def add_leptauh_ABCD_categories(config: od.Config) -> None:
    add_category(config,name="lepA",  id=200,  selection="cat_ss_iso2_bveto_lowmt",      label="A",   tags={"ss","iso2",   "bveto",  "lowmt" })
    add_category(config,name="lepB",  id=400,  selection="cat_ss_noniso2_bveto_lowmt",   label="B",   tags={"ss","noniso2","bveto",  "lowmt" })
    add_category(config,name="lepA0", id=600,  selection="cat_os_iso2_nobveto_lowmt",    label="A0",  tags={"os","iso2",   "nobveto","lowmt" })
    add_category(config,name="lepB0", id=800,  selection="cat_os_noniso2_nobveto_lowmt", label="B0",  tags={"os","noniso2","nobveto","lowmt" })
    add_category(config,name="lepA1", id=1000, selection="cat_os_iso2_bveto_highmt",     label="A1",  tags={"os","iso2",   "bveto",  "highmt"})
    add_category(config,name="lepB1", id=1200, selection="cat_os_noniso2_bveto_highmt",  label="B1",  tags={"os","noniso2","bveto",  "highmt"}) 
    add_category(config,name="lepD",  id=1400, selection="cat_os_iso2_bveto_lowmt",      label="D",   tags={"os","iso2",   "bveto",  "lowmt" })
    add_category(config,name="lepC",  id=1600, selection="cat_os_noniso2_bveto_lowmt",   label="C",   tags={"os","noniso2","bveto",  "lowmt" })

@call_once_on_config()
def add_lept_cp_categories(config: od.Config) -> None:
    add_category(config, name="pi_2",        id=2,  selection="cat_pi_2",      label=r"$\pi$",                 tags={"tau2pi"    }) # h2 -> pi
    add_category(config, name="rho_2",       id=4,  selection="cat_rho_2",     label=r"$\rho$",                tags={"tau2rho"   }) # h2 -> rho
    add_category(config, name="a1dm2_2",     id=6,  selection="cat_a1dm2_2",   label=r"$a_{1}(2p-2\pi^{0})$",  tags={"tau2a1DM2" }) # h2 -> a1
    add_category(config, name="a1dm10_2",    id=8,  selection="cat_a1dm10_2",  label=r"$a_{1}(3p-0\pi^{0})$",  tags={"tau2a1DM10"}) # h2 -> a1
    add_category(config, name="a1dm11_2",    id=10, selection="cat_a1dm11_2",  label=r"$a_{1}(3p-1\pi^{0})$",  tags={"tau2a1DM11"}) # h2 -> a1

    

@call_once_on_config()
def add_etau_mutau_real_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("etau"),
                    config.get_category("mutau")],
        "TorF"   : [config.get_category("real_2")],
        #"abcd"   : [config.get_category("lepA"),  config.get_category("lepB"),
        #            config.get_category("lepA0"), config.get_category("lepB0"),
        #            config.get_category("lepA1"), config.get_category("lepB1"),
        #            config.get_category("lepC"),  config.get_category("lepD")],
        # "abcd"   : [config.get_category("lepC"),  config.get_category("lepD")],
        # "cp"     : [config.get_category("pi_2"),
        #            config.get_category("rho_2"),
        #          config.get_category("a1dm2_2"),
                #    config.get_category("a1dm10_2"),
                #    config.get_category("a1dm11_2")],
    }
    logger.info("etau_mutau_real_categories")
    create_nested_categories(config, categories)

    
@call_once_on_config()
def add_etau_mutau_fake_categories(config: od.Config) -> None:
    categories = {
        "channel": [config.get_category("etau"),
                    config.get_category("mutau")],
        "TorF"   : [config.get_category("fake_2")],
        #"abcd"   : [config.get_category("lepA"),  config.get_category("lepB"),
        #            config.get_category("lepA0"), config.get_category("lepB0"),
        #            config.get_category("lepA1"), config.get_category("lepB1"),
        #            config.get_category("lepC"),  config.get_category("lepD")],
        "abcd"   : [#config.get_category("lepC"),
                    config.get_category("lepD")],
        #"cp"     : [config.get_category("pi_2"),
        #            config.get_category("rho_2"),
        #            config.get_category("a1dm2_2"),
        #            config.get_category("a1dm10_2"),
        #            config.get_category("a1dm11_2")],
    }
    logger.info("etau_mutau_fake_categories")
    create_nested_categories(config, categories)
    

# ########################### #
# Add extra common categories #
# ########################### #
    
@call_once_on_config()
def add_common_categories(config):

    add_category(config,
                 name="incl",
                 id=900000,
                 selection="cat_incl",
                 label=r"$Inclusive$",
                 tags={"incl"})
    add_category(config, name="real_1", id=10000, selection="cat_real_1", label=r"lep1 real", tags={"tau1isRealMC"})
    add_category(config, name="fake_1", id=20000, selection="cat_fake_1", label=r"lep1 fake", tags={"tau1isFakeMC"})
    add_category(config, name="real_2", id=30000, selection="cat_real_2", label=r"lep2 real", tags={"tau2isRealMC"})
    add_category(config, name="fake_2", id=40000, selection="cat_fake_2", label=r"lep2 fake", tags={"tau2isFakeMC"})

    add_category(config, name="has_0j", id=51, selection="cat_0j", label=r"has 0jets",          tags={"has0j"})
    add_category(config, name="has_1j", id=53, selection="cat_1j", label=r"has 1jets",          tags={"has1j"})
    add_category(config, name="has_2j", id=55, selection="cat_2j", label=r"has 2 or more jets", tags={"has2j"})

    add_category(config, name="ipsig_0to1_1",   id=57, selection="cat_ipsig_0to1_1",   label=r"0 < IP_significance < 1", tags={"tau1_ipsig_0to1"})
    add_category(config, name="ipsig_1toany_1", id=59, selection="cat_ipsig_1toany_1", label=r"IP_significance > 1",     tags={"tau1_ipsig_1toany"})

    
# ################### #
# test categorization #
# ################### #
@call_once_on_config()
def add_test_categories(config: od.config) -> None:

    categories = {
        "channel": [config.get_category("tautau")],
        "TorF"   : [config.get_category("real_1")],
        "abcd"   : [
            config.get_category("hadA"),   config.get_category("hadB"),
            config.get_category("hadA0"),  config.get_category("hadB0"),
            config.get_category("hadC0"),  config.get_category("hadD0"),
            config.get_category("hadC"),   config.get_category("hadD"),
        ],
        "IP"     : [config.get_category("ipsig_0to1_1"), config.get_category("ipsig_1toany_1")],
        "cp"     : [config.get_category("pi_rho_1")],
    }

    logger.info("test_categories")
    create_nested_categories(config, categories)
    

    
# ################### #
# main categorization #
# ################### #

@call_once_on_config()    
def add_categories(config: od.Config) -> None:
    """
    Adds all categories to a *config*.
    """
    add_common_categories(config)

    add_category(config,
                 name="etau",
                 id=100000,
                 selection="cat_etau",
                 label=r"$e\tau_{h}$",
                 tags={"etau"})
    add_category(config,
                 name="mutau",
                 id=200000,
                 selection="cat_mutau",
                 label=r"$\mu\tau_{h}$",
                 tags={"mutau"})
    add_category(config,
                 name="tautau",
                 id=400000,
                 selection="cat_tautau",
                 label=r"$\tau_{h} \tau_{h}$",
                 tags={"tautau"})

    # add_leptauh_ABCD_categories(config)
    # add_lept_cp_categories(config)
    add_hadtauh_ABCD_categories(config)
    add_hadr_cp_categories(config)
    
    # add_etau_mutau_real_categories(config)
    #add_etau_mutau_fake_categories(config)
    add_tautau_real_categories(config)
    # add_tautau_fake_categories(config)
    
    # add_test_categories(config)
