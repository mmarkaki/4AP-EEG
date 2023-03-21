# 4AP-EEG 
df/f

The dF/F of individual ROI traces is calculated per segment by subtracting the mean of the lowest 10 percent of the fluorescence values within the segment (i.e., baseline) from each fluorescence value, divided by that baseline:

dF/F=(F-10th Percentile of F)/ 10th Percentile of F                [1]

4-AP-induced activity is characterized by large and relatively long-lasting deviations of the raw fluorescence from baseline (Wenzel M, et al., 2017). To identify a proper baseline avoiding the impact of such events, we computed the interquartile range (IQR) value of the (PMT-corrected) fluorescence in every uninterrupted epoch of acquisition (i.e., segment) in vehicle for each mouse and use their maximum value over all segments as a so-called “IQR-threshold”. The 4-AP segments with IQR value below the IQR threshold for that mouse are identified as quiet. The average 10th percentile of the PMT-corrected fluorescence of the quiet segments is used as the baseline for the less quiet segments in the dF/F estimation below. 

z score

To account for the different noise levels in different ROIs, we generate a “synthetic” distribution (per ROI) by “mirroring” the dF/F values of the frames with dF/F values less than the median of their observed distribution per the entire recording. We then compute the standard deviation of that (synthetic) distribution dist(dF/F) and use that to calculate the z-scored dF/F as follows:

z scored dF/F=(dFF-mean(dFF))/std(dist(dF/F)) 						          [2]

If the z score estimation is performed per segment, then the choice  of a different segment baseline for the df/f estimation in the previous step, has no effect on the z scores. That is, we can simply apply equation [1] for the df/f estimation in every segment, before proceeding to the z score estimation also in every segment.

Local plateaus

We develop a procedure to identify the epochs of significant fluorescence activity in the post 4-AP recording: we calculate noise intervals for each 4AP neuron based on the largest z scores of vehicle neurons. For that, we find the frames with z scores less than the 99.9th percentile of the vehicle neuron's z score. These frames are the noise intervals. We count the consecutive frames of the noise intervals. If there are intervals with less than [a fixed number of] consecutive frames, they are no longer considered as noise interval frames. We count inter-arrival time of remaining noise intervals of the previous step. If inter-arrival time is small (less than 100 frames), then it will no longer be considered as an inter-arrival interval but it will get concatenated together with the prior and following noise interval. Inter-arrivals of larger duration remain the same. Define as (local) valleys the “cleaned up" noise intervals of the previous step. Define as (local) plateaus the in-between frames. Plateaus correspond to epochs of significant activity. Create a time series for each neuron with 0 at frames of valleys and 1 at plateau frames (this is useful for the global plateau/valley definition as well as the matching of global and local plateaus). The concatenation and elimination thresholds have been defined based on the aforementioned duration threshold. 

Global plateaus

A “global plateau” is then defined, starting at the frame when a small percentage of neurons are at their own local plateau (i.e., engagement threshold) and ends at the first frame when the aforementioned neuronal engagement threshold is crossed.  The engagement threshold is selected to be sufficiently small for being able to capture well the entire evolution of the burst of activity, from its start until its end, consistently across mice (see Figures ...). The concatenation and elimination process is also applied on global plateaus, as in local plateaus. The local and global plateaus serve as “auxiliary” markers of the periods of abnormal events in order to characterize the neuronal activity and recruitment process in the ictal events.
After defining the plateaus and valleys on a single neuron level it is time to define them for the entire neuronal population as well. The population or global events are based on the aggregation of the single neuron events: We sum the time series per frame in order to create the aggregate signal based on which we will define the global valleys and plateaus. We find the frames with less than a percentage of neurons having local plateaus simultaneously in order to define the global level noise intervals. We do steps of the procedure followed for single neurons to end up with global valleys and plateaus respectively.
For every global plateau, we choose a “reference window” that has a common start across all these neurons with a local plateau that corresponds to the specific global plateau. This global plateau will be used for assessing the order of recruitment.

Savitzky-Golay (S-G) smoothing

Motion-correction algorithms demonstrate limited effectiveness in correcting tissue pulsation artifacts during high-speed resonant scanning as compared to optical recording involving conventional galvanometers (Wenzel M. et al., 2017). To attenuate this issue and facilitate the fit of local plateaus, the z-scored dF/F of individual neurons was filtered with a Savitzky-Golay (S-G) smoothing envelope that is tolerant of sharp changes in the filtered signal (A. Savitzky and M.J.E. Golay, 1964). The minimal duration of the local plateaus in the vehicle was used to set the window size along with an order 4 local polynomial, which is typically used in S-G smoothing envelopes and works well in different domains. A sensitivity analysis on the impact of 3 window sizes at least equal to the minimal duration of the local plateaus in the vehicle  (101, 211 and 511 frames) yielded the choice of 211 frames as a compromise between smoothing and loss of details during the fit of local plateaus with a sigmoid.

Sigmoid fitting of local plateaus

We proceed to apply a sigmoid fitting of each local plateau, based on the S-G smoothed z-scored dF/F signal across all neurons. The fitting window starts at the first frame of the local plateau of the neuron. If the start of the corresponding segment is nearest, then we take this as the start of the sigmoid fitting. It ends at the maximum amplitude of the detected peak   extended by a small margin, estimated on the smoothed z-scored dF/F signal of the neuron. Specifically, this margin is the frame immediately after the 95th percentile margin -  that follows the detected peak. The start of the sigmoid fitting window is estimated by “mirroring” to the left the frames up to the maximum amplitude with respect to the onset of the local plateau - unless the segment start is earlier. 
The sigmoidal fitting of the local plateaus is defined by three parameters, namely, the maximum amplitude Rmax, the time at which the amplitude is at half maximum S1/2, and a parameter related to the slope at half-maximum amplitude 

rSeff =Rmax/(1+exp(-α(Seff-S1/2))), >0, 

where parameter α is related to the slope at half-maximum amplitude according to the formula 

Slope1/2=αRmax/4

The optimization of parameters of the sigmoid function is performed using the non-linear least square regression nlinfit function of MATLAB 2017a.  The coefficient of determination R2 is used to evaluate the goodness of the fit between the smoothed z-scored dF/F signal of the neuron and the fitted one. Only the fittings with at least R2 equal to 0.8 are taken into account; hence the local plateaus were “filtered” on the basis of good sigmoid temporal fitting. We evaluated mean and standard errors over uniform histogram bins of sigmoid parameters across mice in order to depict their variability. 

We can use the start of the local plateau for the estimation of the order of engagement after we have filtered out the local plateaus with poor sigmoid fitting. We must calibrate local plateaus to the start of the respective population plateau (burst) for defining the order of engagement of each neuron in this specific burst.

Categorizing cells into temporal quartiles and comparing the variance of the distance of each cell with the quartile spatial mean and with the 
variance of the distance to the spatial mean of all cells:
y1=[bigDistB,bigDistG,bigDistY,bigDistR]%distance of each cell to the quartile (blue is the closest, red is the furthest) spatial mean
y2=[bigDistBU,bigDistGU,bigDistYU,bigDistRU]%vs distance of the same cell to the spatial mean of all cells 
y=[y1;y2]
group=[ones(1,105),2*ones(1,231),3*ones(1,224),4*ones(1,249)];%numbers refer to the number of cells belonging to each spatial quartile in all mice
[P,ANOVATAB,STATS] = anova1(y,group)
[c,m,h,nms] = multcompare(STATS,'display','off');


EEG CWT and PSD estimation
