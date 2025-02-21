# nmr-cdiff
Supporting scripts and data for the analysis of real-time metabolism using <sup>13</sup>C NMR spectra of live _C. difficile_ cell cultures, [as published in Nature Chemical Biology (https://doi.org/10.1038/s41589-023-01275-9).](https://doi.org/10.1038/s41589-023-01275-9). Updated from original work to be automated, modular, and accomadate processing of N15 data.

## Folder structure:
> [data](data) - datasets, including the updated `icdf843` metbolic model and the dFBA parameters  
> [etc](etc) - contains requirements.txt for a virtual environment to run the python code  
> [scripts](scripts) - all scripts used in analysis  
> > [process](scripts/process) - a python script and supporting csh scripts to process the raw NMR files
> > [dfba.py]() - a python script supporting dynamic FBA modeling constrained by NMR data
 
## Installation
Clone this repository with `git clone https://github.com/Massachusetts-Host-Microbiome-Center/nmr-cdiff.git`

### Dependencies
1. The remaining analyses require python 3.8+ and a collection of packages. We have provided a `pip freeze` output [requirements.txt](etc/requirements.txt) that includes all of the python dependencies. The simplest way to load these dependencies is by creating a python virtual environment with `venv`, then installing all of the dependencies with `python -m pip install -r requirements.txt`.
2. [NMRpipe](https://www.ibbr.umd.edu/nmrpipe/install.html) is required to run the NMR processing scripts in process.py and any other functions that call them. Please follow the installation instructions for NMRpipe carefully.

## Tutorial
Follow the READMEs laid out in order
1. [1.installing_nmr_processing_environment.md](1.installing_nmr_processing_environment.md)
2. [2.preprocessing_data.md](2.preprocessing_data.md)
3. [3.dfba.md](3.dfba.md)
4. [4.compare_kinetics.md](4.compare_kinetics.md)

## Building Bifermentans metabolic model
5.building_pbi_model.md

## License
This distribution is available under the [Apache License, Version 2.0](LICENSE).
