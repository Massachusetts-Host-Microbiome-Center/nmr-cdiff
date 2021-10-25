# Plot NMR data using MATLAB
## Functions
1. getScaledCoeffs – calculate the time-scaled logistic coefficients to be used dFBA analyses.
2. getTimescale – calculate the timing and rate of metabolism using trajectory of reference peak.
3. plotAlaTraces –  plots NMR traces of the alpha carbon alanine peaks, with splitting lines (Fig. 4A,C).
4. plotAllScaled – scales the time axis of the stack to align with the onset time and rate of metabolism, then plots 13C stacks with colored peaks, and 1D endpoint spectra (Fig. 1A,C,E; Fig. S1).
5. plotCombinedPanel - plot a logistic fit of all metabolites in a single condition (Fig. 1B,D,F).
6. plotFluxBidirectional - stacked area plot, showing proportions of flux for a metabolite. Takes data from "ac_prop" (see dFBA) (Fig. 2B-E).
7. plotGroupedReactions - Plot line plot with error regions, for select dFBA results (Fig. 2A).
8. plotRegionsAla – plot 13C waterfall with colored peaks at Ala-a13C (Fig. 4B,D).
9. plotRegionsLeu – like plotAllScaled, with special considerations for leucine.

## Usage
Some of these functions, especially plotAllScaled, use a large amount of memory. For best results, call the functions in batch mode from the command line.
