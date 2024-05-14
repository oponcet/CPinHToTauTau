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

### A quick guide to create and link properly your DesyTau/CPinHToTauTau working repository 

1. Go to your github and fork the [DesyTau/CPinHToTauTau](https://github.com/DesyTau/CPinHToTauTau) and the [uhh-cms/cmsdb](https://github.com/uhh-cms/cmsdb) repositories.

2. From your workspace do:
```bash
git clone https://github.com/your_github_username/CPinHToTauTau.git
```
3. Inside CPinHTauTau repository do:
```bash
git submodule update --init --recursive
```
4. Then set cmsdb as follow:
```bash
cd CPinHToTauTau/modules/cmsdb
git remote remove origin
git remote add desytau https://github.com/DesyTau/cmsdb.git (common repository)
git remote add upstream https://github.com/uhh-cms/cmsdb.git (original repository)
git remote add origin git@github.com:your_github_username/cmsdb.git (personal repository, the link here is the ssh one)
```
5. Check if everiything is in place by:
```bash
git remote -v
```
You should see something like this:
```bash
desytau https://github.com/DesyTau/cmsdb.git (fetch)
desytau https://github.com/DesyTau/cmsdb.git (push)
origin  git@github.com:jmalvaso/cmsdb.git (fetch)
origin  git@github.com:jmalvaso/cmsdb.git (push)
upstream        https://github.com/uhh-cms/cmsdb.git (fetch)
upstream        https://github.com/uhh-cms/cmsdb.git (push)
```
6. Besides, you need to pull the current cmsdb version from DesyTau/cmsdb:
```bash
git pull desytau master
```
7. Everything should be in a good shape now, go back to your local CPinHToTauTau repository:
 ```bash
cd ../..
```
and run the setup.sh:
```bash
source setup.sh dev
```
Enjoy your perfectly github-linked local repository :)
