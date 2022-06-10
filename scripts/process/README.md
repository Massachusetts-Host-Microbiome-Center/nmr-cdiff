# NMR pre-processing]
## Note
See [README.md](../README.md) in the base directory for dependencies before attempting to run this script. Be sure to run process.py or any other scripts that call its functions in a c-shell with NMRpipe installed.
## Usage
To process a run:
1. Place a text file `cfg_1H.txt` or `cfg_13C.txt` in the raw NMR file directory, containing the reference ppm shifts and names of the expected chemicals, one per line, tab-separated. See below for an example.  
2. Ensure [calibrations.txt](calibrations.txt) contains a row for your run.  
3. Navigate to the directory containing the raw NMR files from bruker using `cd`.  
4. Process the run using `python /path/to/process.py [isotope] [init]`.  
    - `[isotope]`: the isotope label (13C or 1H). Default=13C.  
    - `[init]`: the ID of the initial spectrum, used as the time anchor. Default=11.  
5. The script will prompt you at several points. Use the nmrglue phase corrector on each spectrum to adjust as necessary. If prompted to assign peaks to compounds, make note of the peak contributions and close the plot window. Then type the indices of each peak belonging to a compound as a whitespace-separated list.  

Example of `cfg_13C.txt`:
```
177.336 Proline
64.0059	Proline
48.8949	Proline
31.632	Proline
26.426	Proline
185.746	5-aminovalerate
42.1713	5-aminovalerate
39.4304	5-aminovalerate
29.2843	5-aminovalerate
25.2451	5-aminovalerate
```
[HMDB](https://hmdb.ca/) is a good source for reference chemical shifts.

To run dFBA analyses:
1. After all runs are processed, evaluate the conditions for the best scaling approach. If necessary, add the substrate to `substrate_map` and add a condition to define `sfacts` for the substrate in `trajectories.py`. If an NMR-standards approach is to be taken (as we did for 13C-glucose), process the standards like an NMR run and provide "Concentrations.xlsx" giving the input concentrations of each standard in the same shape as the "areas" tab of the output from `process.py`.
2. Run the dFBA analyses using `python /path/to/dfba.py [method] [/path/to/run1] [/path/to/run2] ... `  
    - `[method]` is fba, pfba, loopless, or fva
