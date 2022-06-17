# Tutorial for NMR processing pipeline
This tutorial assumes a working knowledge of the command line. The instructions were written for Mac users but could be adapted for Windows systems with small modifications. Each section includes an estimate of how long it will take.  

1. [Preparing to run the pipeline (20 minutes)](TUTORIAL.md#preparing-to-run-the-pipeline-20-minutes)
2. [Processing the NMR runs (15 minutes per run)](TUTORIAL.md#processing-the-nmr-runs-15-minutes-per-run)
3. [Running the dFBA analyses (5 minutes)](TUTORIAL.md#running-the-dfba-analyses-5-minutes)

## Preparing to run the pipeline (20 minutes)

### Installing the dependencies

Make sure all requirements are satisfied in the "Installation" and "Dependencies" sections of [README.md](README.md). A c-shell environment (csh/tcsh) must be activated when installing NMRPipe and working with `process.py`: 
1. Install [NMRPipe](https://www.ibbr.umd.edu/nmrpipe/install.html). On a Mac, the installation process should look something like this:
    ```
    csh
    cd
    mkdir nmr
    cd nmr

    curl -O https://www.ibbr.umd.edu/nmrpipe/install.com
    curl -O https://www.ibbr.umd.edu/nmrpipe/binval.com
    curl -O https://www.ibbr.umd.edu/nmrpipe/NMRPipeX.tZ
    curl -O https://www.ibbr.umd.edu/nmrpipe/s.tZ
    curl -O https://www.ibbr.umd.edu/nmrpipe/dyn.tZ

    chmod a+r  *.tZ *.Z *.tar
    chmod a+rx *.com

    ./install.com
    ```
    Check your `.cshrc` file for the initialize commands.
    ```
    cd
    cat .cshrc
    ```
      The file should contain the following lines:
    ```
    if (-e /path/to/nmr/com/nmrInit.mac11_64.com) then
       source /path/to/nmr/com/nmrInit.mac11_64.com
    endif
    ```
    If the lines are missing or the .cshrc file does not exist, open the file `README_NMRPIPE_USERS` in the `nmr` folder you just created and find the lines that match the ones above. Copy them into .cshrc file. Then log out and log back in again.  

2. Download the repository with `git clone https://github.com/Massachusetts-Host-Microbiome-Center/nmr-cdiff.git` or using the GitHub Desktop app.  
3. Set up the virtual environment. From the `nmr-cdiff` base directory:
    ```
    cd etc
    python3 -m venv env
    source env/bin/activate
    python3 -m pip install -r requirements.txt
    ```

### Test data

One short unprocessed NMR run and two processed runs are included in the [test data](data/test) folder:
 - [20210519_13CGlc](data/test/20210519_13CGlc) : a short unprocessed glucose run  
 - [20210322_13CLeu](data/test/20210322_13CLeu) : a short processed leucine run  
 - [20220421_13CGlc_standards](data/test/20220421_13CGlc_standards) : a processed dataset of standard solutions supporting the dFBA analyses

We will process the glucose run and simulate dFBA. For convenience, we have included `cfg_13C.txt` and `cfg_1H.txt` in the run directory, which contain reference chemical shifts for the substrate and expected products. For each substrate, these files must be prepared before `process.py` is run.

## Processing the NMR runs (15 minutes per run)
### Execute process.py
We will begin by processing the NMR datasets using the [process.py](scripts/process/process.py) script. Navigate to the run directory and run the processing script.  
```
csh
source nmr-cdiff/etc/env/bin/activate.csh 
cd nmr-cdiff/data/test/20210519_13CGlc
python ../../../scripts/process/process.py 13C 52
```
The program will identify the 13C scans and use scan 52 as the timestamp anchor.

### Find the calibration parameters
Let's first find some calibration parameters to make processing go faster. 

1. You will first be prompted for phase correction.  
    <img width="752" alt="Phasing_unphased" src="https://user-images.githubusercontent.com/29278926/174334460-d3786f90-8032-415b-8d8c-f78a4e68eea2.png">  
    Slide the p0 and p1 bars until the spectrum is well-phased.  

    <img width="752" alt="Phasing_phased" src="https://user-images.githubusercontent.com/29278926/174334459-bed1879d-a99d-4533-90dc-30688a1cb535.png">  

    Press the "Set Phases" button, and close the window. For this spectrum, **p0=70.3** and **p1=99.9** were appropriate values.  

2. Next, you will be prompted to calibrate the reference shift.
    - Find a well-separated peak in the reference 13C spectrum of Glucose. The peak at **72.405 ppm** appears suitable in [this reference spectrum from HMDB](https://hmdb.ca/spectra/nmr_one_d/166522).
    - In the plot that pops up from process.py, find the corresponding peak in the spectrum, which may be slightly shifted from 72.405 ppm. Use the zoom tool to locate the peak and cursor to find the chemical shift at the center of the peak. Write down the chemical shift and close the window.  

        <img width="752" alt="Calib1" src="https://user-images.githubusercontent.com/29278926/174334449-7f676cb0-5b16-4b1f-a443-c6777a0fb844.png">  

    - In the command line, type the chemical shift that you just wrote down. Hit return. We found that **69.95 ppm** was the approximate experimental shift for this peak.

        <img width="752" alt="Calib2" src="https://user-images.githubusercontent.com/29278926/174335367-003a83b1-73f0-4e8d-bed6-ec5ff97dde8d.png">

These parameters will be stored to aid with processing each remaining spectrum in the run.

### Process the 13C spectra
Processing will now begin, and should take several minutes. As each spectrum is processed, you will be prompted for input one or more times.  

1. The phase correction window will appear. Adjust the phase correction as before if necessary and close the window.  
3. At times, manual peak assignment may be required if peaks are close to each other. If this is the case, a plot will show where each peak is labeled by a numbered index.In the terminal window, the conflicting reference shifts will be listed.  

    <img width="643" alt="Splitting1" src="https://user-images.githubusercontent.com/29278926/174334468-1e93f9aa-3a62-4e6e-ad95-8811a0c9f4a1.png">  

    Make note of which peaks belong to which compounds and close the window. Then, in the command line window, list the indices belonging to each compound separated by whitespace.  
    
    <img width="752" alt="Splitting2" src="https://user-images.githubusercontent.com/29278926/174335504-d02fb9b9-7108-4635-a8da-23e6b7b2d5dc.png">  

6. The SciPy routines executing the peak-fitting may occasionally print warnings. These can be ignored.  
7. After peak picking and fitting is completed for each spectrum, the spectrum will be plotted with the curve-fit and peak assignments overlaid. Advanced users may wish to inspect this spectrum for proper curve-fitting. If curve-fitting is consistently poor, the parameters `r` and `sep` may be adjusted in the calls to `peak_fit` and `prominent` in `process.py`.  

    <img width="752" alt="Proc1" src="https://user-images.githubusercontent.com/29278926/174334464-0e00b470-0a18-4633-b916-d4572986e90b.png">  

When processing is finished, the spectrum stacks will be plotted and trajectories the peaks in `cfg_13C.txt` will be displayed. Next, the processed output will be written to `20210519_13CGlc_13C.xlsx`.

<img width="752" alt="Proc2" src="https://user-images.githubusercontent.com/29278926/174334465-a8c58721-27aa-4bd5-872a-f06f8aa73d1b.png">

### Process the 1H spectra
Now, run the processing script for the 1H spectra.
```
python ../../../scripts/process/process.py 1H 52
```
Follow the prompts as before. We used p0=155.8, p1=50.2, reference shift=1.9059, experimental shift=1.9046.

<img width="752" alt="Proc3" src="https://user-images.githubusercontent.com/29278926/174334467-2459bbe8-c640-4bff-86c6-a4cff6e99696.png">

### Plotting the runs
A fully processed run can be plotted as a stack, as seen in Fig. 2a,c,e, using the MATLAB function [`plotStacks.m`](scripts/MATLAB/plotStack.m). This can be accomplished with relative ease using the script [`plot_stack.m`](scripts/MATLAB/call/plot_stack.m). Simply add a case for the run and set the variable `met`. You may also need to add a working filepath using the funciton `addpath`. The script is best run in batch mode from the command line:
```
matlab -nodisplay -nodesktop -batch "run('plot_stack.m');exit;"
```

### Additional processed runs
- The run [`20210322_13CLeu`](data/test/20210322_13CLeu) has already been processed for convenience.  
- The run [`20220421_13CGlc_standards`](data/test/20220421_13CGlc_standards) is a processed dataset of standard solutions supporting concentration estimates in the glucose run. The file `Concentrations.xlsx` in this directory contains the concentration of each compound in the processed spectra.

## Running the dFBA analyses (5 minutes)
### Overview of the script
Now that the runs are processed, we will simulate dFBA. The dFBA script will perform the following functions:  
1. Synchronize the NMR runs using the isocaproate trajectory from the 1H spectra.
2. Estimate concentration curves of the substrates and products using standard solutions.
    - The standards for glucose are included in the run `20220421_13CGlc_standards` and will be used to create standard curves that estimate relative NMR snesitivity of each compound with respect to glucose.
    - The expected ratios of isocaproate and isovalerate to leucine were estimated from gas chromatography experiments, and are included directly in  [`trajectories.py`](scripts/process/trajectories.py).
    - Plots of these curves as seen in Fig. 2b,d,f will be saved in each run folder.
3. Simulate dFBA using the estimated concentration curves as constraints for exchange fluxes.  

### Run the dFBA analyses
Run the dFBA script, [`dfba.py`](scripts/process/dfba.py):
```
cd ~/nmr-cdiff
python scripts/process/dfba.py fba data/test/20210519_13CGlc_13C data/test/20210322_13CLeu_13C
```
The first argument, `fba`, tells the script that we are computing a standard FBA solution at each time point. Next, we specify the directories containing the runs.

By default, we will only track a defined set of reactions and metabolites. These tracked elements can be modified by editing the `tracked_reactions` and `tracked_metabolites` variables near the top of `dfba.py`.

### dFBA output
The dFBA program will display some text output and graphs with the simulation results. Output from tracked reactions and metabolites will be written to `data/fluxes.xlsx`, `data/met_fluxes.xlsx`, and `data/<met>flux_in.txt` and `data/<met>flux_out.txt`. These files support MATLAB plotting functions.

The logistic concentration curves should appear as follows:

<img width="362" alt="dfba1" src="https://user-images.githubusercontent.com/29278926/174334454-8b1242f2-aacd-4e83-b634-4b191be7154d.png">
<img width="1552" alt="dfba2" src="https://user-images.githubusercontent.com/29278926/174334456-e5621b9c-ac7d-4ecc-a965-21d9ed6ea5dc.png">

Our dFBA output was:

<img width="1552" alt="dfba3" src="https://user-images.githubusercontent.com/29278926/174334458-77a59d69-68e0-41fb-bd0a-94de486ebf48.png">

