# Tutorial for NMR processing pipeline
This tutorial assumes a working knowledge of the command line. The instructions were written for Mac users but could be adapted for Windows systems with small modifications.  

## Preparing to run the pipeline (15 minutes)

### Installing the dependencies

Make sure all requirements are satisfied in the "Installation" and "Dependencies" sections of [README.md](README.md): 
1. Install [NMRPipe](https://www.ibbr.umd.edu/nmrpipe/install.html). On a Mac, the installation process should look something like this:
    ```
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

Two short NMR runs are included in the [test data](data/test) folder:
 - 20210103_13CGlc : a short glucose run  
 - 20210322_13CLeu : a short leucine run  

We will process these runs and run a dFBA simulation.

## Processing the NMR runs (30 minutes)
We will begin by processing the NMR datasets using the process.py script. Let's first find some calibration parameters to make processing go faster. Navigate to the run directory and run the processing script.  
```
csh
source nmr-cdiff/etc/env/bin/activate.csh 
cd nmr-cdiff/data/test/20210103_13CGlc
python ../../../scripts/process/process.py 13C 11
```
The program will identify the 13C scans and use scan 11 as the timestamp anchor.  

You will first be prompted for phase correction. Slide the p0 and p1 bars until the spectrum is well-phased. Press the "Set Phases" button, and close the window. For this spectrum, p0=84 and p1=81 were ideal values.  

Next, you will be prompted to calibrate the reference shift.
- Find a well-separated peak in the reference 13C spectrum of Glucose. The peak at **72.405 ppm** appears suitable in [this reference spectrum from HMDB](https://hmdb.ca/spectra/nmr_one_d/166522).
- In the plot that pops up from process.py, find the corresponding peak in the spectrum, which may be slightly shifted from 72.405 ppm. Use the zoom tool to locate the peak and cursor to find the chemical shift at the center of the peak. Write down the chemical shift and close the window.
- In the command line, type the chemical shift that you just wrote down. Hit return.

These parameters will be stored to aid with processing each remaining spectrum in the run.

Processing will now begin, and should take several minutes. As each spectrum is processed, the phase correction window will appear. Adjust the phase correction as before if necessary and close the window.

At times, manual peak assignment may be required if peaks are close to each other. If this is the case, a plot will show where each peak is labeled by a numbered index. In the terminal window, the conflicting reference shifts will be listed. Make note of which peaks belong to which compounds and close the window. Then, in the command line window, list the indices belonging to each compound separated by whitespace.  
-- TO BE CONTINUED --
