# coding: utf-8

"""
Configuration of the CPinHToTauTau analysis.
"""
import law
import order as od

# ------------------------ #
# The main analysis object #
# ------------------------ #
analysis_httcp = ana = od.Analysis(
    name="analysis_httcp",
    id=1,
)

# analysis-global versions
# (see cfg.x.versions below for more info)
ana.x.versions = {}

# files of bash sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.bash_sandboxes = ["$CF_BASE/sandboxes/cf.sh"]
default_sandbox = law.Sandbox.new(law.config.get("analysis", "default_columnar_sandbox"))
if default_sandbox.sandbox_type == "bash" and default_sandbox.name not in ana.x.bash_sandboxes:
    ana.x.bash_sandboxes.append(default_sandbox.name)

# files of cmssw sandboxes that might be required by remote tasks
# (used in cf.HTCondorWorkflow)
ana.x.cmssw_sandboxes = [
    "$CF_BASE/sandboxes/cmssw_default.sh",
]

# config groups for conveniently looping over certain configs
# (used in wrapper_factory)
ana.x.config_groups = {}

# ------------- #
# setup configs #
# ------------- #

# ------------------------ Run2 UL 2017 ----------------------- #
from httcp.config.configs_run2ul_SR import add_config as add_config_run2ul_SR
from cmsdb.campaigns.run2_2017_nano_local_v10 import campaign_run2_2017_nano_local_v10
add_config_run2ul_SR(
    analysis_httcp,
    campaign_run2_2017_nano_local_v10.copy(),
    config_name=f"{campaign_run2_2017_nano_local_v10.name}_GS",
    config_id=1,
)
add_config_run2ul_SR(
    analysis_httcp,
    campaign_run2_2017_nano_local_v10.copy(),
    config_name=f"{campaign_run2_2017_nano_local_v10.name}_limited",
    config_id=2,
    limit_dataset_files=1,
)
# ------------------------------------------------------------- #

# ------------------------ Run2 UL 2017 ----------------------- #
from httcp.config.run2_UL2017 import add_run2_UL2017
from cmsdb.campaigns.run2_UL2017_nano_tau_v10 import campaign_run2_UL2017_nano_tau_v10
add_run2_UL2017(
    analysis_httcp,
    campaign_run2_UL2017_nano_tau_v10.copy(),
    config_name=f"{campaign_run2_UL2017_nano_tau_v10.name}_limited",
    config_id=3,
    limit_dataset_files=1)
# ------------------------------------------------------------- #

# -------------------------- Run3 2022 ------------------------ #
from httcp.config.run3_preEE import add_run3_preEE
from cmsdb.campaigns.run3_2022_preEE_nano_tau_v12 import campaign_run3_2022_preEE_nano_tau_v12
add_run3_preEE(analysis_httcp,
                  campaign_run3_2022_preEE_nano_tau_v12.copy(),
                  config_name=f"{campaign_run3_2022_preEE_nano_tau_v12.name}",
                  config_id = 4)
add_run3_preEE(analysis_httcp,
                  campaign_run3_2022_preEE_nano_tau_v12.copy(),
                  config_name=f"{campaign_run3_2022_preEE_nano_tau_v12.name}_limited",
                  config_id = 5,
                  limit_dataset_files=1)
# ------------------------------------------------------------- #


# -------------------------- Run3 2022 ------------------------ #
# run3_2022_postEE_nano_tau_v12_limited
from httcp.config.run3_postEE import add_run3_postEE
from cmsdb.campaigns.run3_2022_postEE_nano_tau_v12 import campaign_run3_2022_postEE_nano_tau_v12
add_run3_postEE(analysis_httcp,
                  campaign_run3_2022_postEE_nano_tau_v12.copy(),
                  config_name=f"{campaign_run3_2022_postEE_nano_tau_v12.name}",
                  config_id = 6)
add_run3_postEE(analysis_httcp,
                  campaign_run3_2022_postEE_nano_tau_v12.copy(),
                  config_name=f"{campaign_run3_2022_postEE_nano_tau_v12.name}_limited",
                  config_id = 7,
                  limit_dataset_files=1)
# ------------------------------------------------------------- #
