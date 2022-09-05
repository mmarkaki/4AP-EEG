function [plateau_start, plateau_end, timeseries_gl] = global_plateaus_init(forConcat,k,timeseriesl)
%for each one of the neurons
% find frames that are in noise interval;
% construct the signal with 0 at noise interval frames and >0 values in
% the rest of the frames
    [cellnum,number_of_frames]=size(timeseriesl);
    timeseries_gl=ones(number_of_frames,1); % here we will keep the 0-1 timeseries
    for_global = sum(timeseriesl,2); 
    
    plateau_start=[];
    plateau_end=[];
    valley_start=[];
    valley_end=[];
    valley_start_c=[];
    valley_end_c=[];
    
    temp = find(for_global<=k*cellnum/100);
    valley_start = temp(1);
    for j=2:length(temp)
        if temp(j)>temp(j-1)+1 & j~=length(temp)
            valley_end = [valley_end; temp(j-1)];
            valley_start = [valley_start; temp(j)];
        elseif j==length(temp)
            valley_end = [valley_end; temp(j)];
        end
    end
    % drop valleys of small duration
    valley_dur = [];
    valley_dur = valley_end - valley_start;
    temp = find(valley_dur<=forConcat);%forConcat 0 leads to more, smaller plateaus
    valley_start(temp)=[];
    valley_end(temp)=[];
    valley_start_c=valley_start(1);
    valley_end_c=valley_end(1);
    counter = 1;
    % if the interval between 2 consecutive valleys is really small
    % discard;small plateaus are omitted this way
    for i=2:length(valley_start)
        temp = valley_start(i)-valley_end(i-1)-2; % valley_start(i)-1 - (valley_end(i-1)+1)
        if temp<=forConcat
            valley_end_c(counter) = valley_end(i);
        else
            valley_start_c=[valley_start_c; valley_start(i)];
            valley_end_c=[valley_end_c; valley_end(i)];
            counter=counter+1;
        end
    end
    % define plateaus
    if valley_start_c(1)>1 % if 1st frame is not valley....
        plateau_start=[plateau_start; 1];
    else
        plateau_start=[plateau_start; valley_end_c(1)+1];
    end
    for i=1:length(valley_start_c) % define intermediate plateaus
        plateau_end=[plateau_end; valley_start_c(i)-1];
        plateau_start=[plateau_start; valley_end_c(i)+1];
    end
    if valley_end_c(end)==number_of_frames % define the end of the last plateau
        plateau_end=[plateau_end; valley_start_c(end)+1];
    else
        plateau_end=[plateau_end; number_of_frames];
    end   
    % timeseries with 0 at valleys and 1 at plateaus/ usefull to visualize
    % plateaus and valleys
    for i=1:length(valley_start_c)
        timeseries_gl(valley_start_c(i):valley_end_c(i))=0;
    end
    plateau_start(1)=[];
    plateau_start(end)=[];
    plateau_end(1)=[];
    plateau_end(end)=[]; 
end

