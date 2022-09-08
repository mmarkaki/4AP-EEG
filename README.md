# 4AP-EEG 
The dF/F of individual ROI traces is calculated then per segment by subtracting the mean of the lowest 10 percent of the fluorescence values within the segment (i.e., baseline) from each fluorescence value, divided by that baseline:

dF/F=(F-10th Percentile of F)/ 10th Percentile of F            

4-AP-induced activity is characterized by large and relatively long-lasting deviations of the raw fluorescence from baseline (Wenzel M, et al., 2017). To identify a proper baseline avoiding the impact of such events, we computed the interquartile range (IQR) value of the PMT-corrected fluorescence in every uninterrupted epoch of acquisition (i.e., segment) in vehicle for each mouse and use their maximum value over all segments as a so-called “IQR-threshold”. The 4-AP segments with IQR value below the IQR threshold for that mouse are identified as quiet. The average 10th percentile of the PMT-corrected fluorescence of the quiet segments is used as the baseline for the less quiet segments in the dF/F estimation below. 
