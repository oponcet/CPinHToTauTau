# CPinHToTauTau Analysis

A Columnflow based analysis framework from IPHC and DESY

<!-- marker-before-logo -->

<div style="text-align: center;">
    <img src="assets/logo.png" alt="Logo" style="width: 400px; height: 220px; display: block; margin: 0 auto;">
</div>

<!-- marker-after-logo -->

### Resources

- [columnflow](https://github.com/columnflow/columnflow/tree/master)
- [law](https://github.com/riga/law)
- [order](https://github.com/riga/order)
- [luigi](https://github.com/spotify/luigi)


### Installation

```sh
git clone --recursive https://github.com/gsaha009/CPinHToTauTau.git
cd CPinHToTauTau

# set environment
source setup.sh hcp
```

Then follow the details:

```sh
CERN username (CF_CERN_USER, default gsaha):  
Local data directory (CF_DATA, default ./data):  
Local directory for installing software (CF_SOFTWARE_BASE, default $CF_DATA/software):  /eos/user/g/gsaha/CPinHToTauTauData                            
Local directory for storing job files (CF_JOB_BASE, default $CF_DATA/jobs):             
Relative path used in store paths (see next queries) (CF_STORE_NAME, default cf_store):  
Default local output store (CF_STORE_LOCAL, default $CF_DATA/$CF_STORE_NAME):  
Local directory for caching remote files (CF_WLCG_CACHE_ROOT, default ''):  
Automatically update virtual envs if needed (CF_VENV_SETUP_MODE_UPDATE, default False):  
Use a local scheduler for law tasks (CF_LOCAL_SCHEDULER, default True):  
Flavor of the columnflow setup ('', 'cms') (CF_FLAVOR, default cms):  
storage element for crab specific job outputs (e.g. T2_DE_DESY) (CF_CRAB_STORAGE_ELEMENT, default ''):  
base directory on storage element for crab specific job outputs (CF_CRAB_BASE_DIRECTORY, default /store/user/$CF_CERN_USER/cf_crab_outputs):
```

This will install the `softwares` in `eos`, but the `data` and `jobs` directories are in `afs`.
In the next step, the big output files have to be stored in `eos`.
So, `law.cfg` needs to be set like here in this repo.

For any computing system, please make sure that this hard-coded (stupid) line is changed to the proper one:

`thisdir` in the config e.g. here: https://github.com/gsaha009/CPinHToTauTau/blob/main/httcp/config/config_run3.py#L33

### Best Practices

 - If you want to check whether your changes are working properly, use the config with the "_limited" tag in the end of your actual config. It will take one ROOT file as input to make the execution faster.
 - For production, if you need to send jobs to any batch system configured in `law`, it is better to run `cd.ReduceEventsWrapper` first. Then it is safe to use `cf.ProduceColumns`. By the end of this task, all possible event and object level corrections have already been applied. So, now usually there are three targets: Plotting, Saving the skimmed and corrected events array (porbably in flat ROOT ntuple) or to produce DataCard. In between, ML training or evaluation task can be used if needed. It is better to follow this `columnflow` [task-graph](https://github.com/columnflow/columnflow/wiki#default-task-graph), always.
 - Now, you may have a ton of datasets and processes. This [script](https://github.com/gsaha009/CPinHToTauTau/blob/main/cf_run.py) might help you to get the commands, which you can copy and paste in a `tmux/screen` session.
   - To start with, you can modify this [`yaml`](https://github.com/gsaha009/CPinHToTauTau/blob/main/yamls/2022PreEE_full.yml) depending on your need. It is recommended to make separate `yaml` for different eras.
   - Then you can run the `cf_run` script like `python cf_run.py -i <yaml/2022PreEE_full.yml or other> -f <ReduceEvents or other>`
   - Notice the version name. If it is dummy, the script will produce some version on it's own, otherwise of you want to use some older version, just specify the name there.
   - And, at the end, you will get the command you want to run.
 - While using `cf.PlotVariables1/2D`, better not to mention all the categories in the command, because it will take ages as it creates the categories in the runtime. So, it would be better to produce plots in several steps. 