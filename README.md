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