# Plot NMR data using MATLAB
## Functions
1. plotCombinedPanel – plot average integrated signal for all compounds, using assembled xlsx file as input. See Figure 1B.
2. plotCompare - plot average integrated signal, comparing two conditions. See Figure 2C.
3. plotFluxArea - plot fluxes of competing reactions, stacked. Takes output from dFBA (`pyr_prop`). See Figure 3C.
4. plotRegions - plot NMR stack as a surface plot with color-labeled ridges. See Figures 1A, 2A, and 2B.

## Usage
Some of these functions, especially plotRegions, use a large amount of memory. For best results, call the functions in batch mode from the command line.
