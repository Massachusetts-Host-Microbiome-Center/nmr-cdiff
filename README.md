# nmr-cdiff
Supporting scripts and data for the analysis of real-time metabolism using <sup>13</sup>C NMR spectra of live _C. diff_ cell cultures.

## Folder structure:
> [data](data) - datasets, including the updated `icdf834` metbolic model and the dFBA parameters  
> [scripts](scripts) - all scripts used in analysis 
> > [dfba](scripts/dfba) - a jupyter notebook with the dFBA implementation  
> > [plot](scripts/plot) - MATLAB scripts for visualizing the NMR data in various plot types  
> > [process](scripts/process) - a python script and supporting csh scripts to process the raw NMR files  
 
## Installation
Clone this repository with `git clone https://github.com/Massachusetts-Host-Microbiome-Center/nmr-cdiff.git`

## Dependencies
The different components of this package have varying dependencies:
1. The dFBA analysis requires an installation of python with [Jupyter](https://jupyter.org/install), [COBRApy](https://opencobra.github.io/cobrapy/), numpy, pandas, and matplotlib. We recommend installing these in a virtual environment using `venv` or `conda`. The notebook was tested and developed with python 3.8.5.
2. The MATLAB functions require an installation of MATLAB. They were developed and tested with MATLAB R2019b.
3. The processing scripts for the raw NMR files require [NMRpipe](https://www.ibbr.umd.edu/nmrpipe/install.html), python 3.8+, and a collection of python packages. We have provided a `pip freeze` output [requirements.txt](etc/requirements.txt) that includes all of the python dependencies. The simplest way to load these dependencies is by creating a python virtual environment with `venv`, then installing all of the dependencies with `python -m pip install -r requirements.txt`.
