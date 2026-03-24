# SpinAlignment Analysis

This directory contains code and resources for the analysis of Spin Alignment.

## Structure

- `steeringFile/`: Contains steering files for running analyses:
    - Belle II: `steeringFile/b2bii_test.py`
    - Belle I: `steeringFile/belle1_steeringFile/SpinAlignment.cc`
    - **Tracking steering files with git tags:**  
        To associate a steering file with a specific root file, use git as follows:
        1. Stage changes: `git add steeringFile`
        2. Commit: `git commit -m "<commit message>"`
        3. Tag the commit: `git tag -a root_produ_<rootFileName> -m "<tag message>"`

- `offline/` : Contains script for offline analysis (plotting, ...): `offline/AnaTest_SpinAlignment.py*`

- `rootFiles/`: Store  symbolic links to input/output ROOT files (avoid duplicating large files). 
When a steering file is needed, create a symlink to it in coordinate directory 
and symlink the produced ROOT output here after the run. Examples:

## Other
- Additional packages used in this analysis are custom-developed and installed via pip. The source directory for these packages is: `https://gitlab.desy.de/wangz731/ANA_TOOLS`.


## Contact

My email : `wang731@mail.ustc.edu.cn`
