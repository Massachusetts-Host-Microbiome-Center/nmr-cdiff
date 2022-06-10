# Plot NMR data using MATLAB
## Functions
1. getColor - retreive colors for visualizations
2. getTimescale - normalize timescale of experiment using 1H spectrum, analogous to synchrnoizers.py
3. load13CData - load Excel output from process.py for plotting
4. plotAlaTraces - plot the splitting patterns of alanine, as in Fig. 5
5. plotFluxBidirectional - plot fluxes of competing reactions, stacked. Takes output from dfba.py (`fluxes.xlsx`). See Figure 3b.
6. plotGroupedRxns - plot fluxes of key reactions, as in Fig. 3a,c-i
7. plotRegionAla - plot NMR stack of alanine time series with color-coded peaks
8. plotStacks - plot NMR stack as a surface plot with color-labeled ridges. See Figure 2a,c,e.

## Usage
Some of these functions, especially plotStacks, use a large amount of memory. For best results, call the functions in batch mode from the command line.
