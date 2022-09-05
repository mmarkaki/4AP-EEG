function [neuron_a] = local_plateaus_neuron(y,for_concat)
%for each one of the neurons
% find frames that are in noise interval;
% construct the signal with 0 at noise interval frames and >0 values in
% the rest of the frames
%y = psm_fl1_yan_no_filt_20(:,n).*(df_f(:,n)); 
    df_f_for_plateaus=y;
    neuron_a = struct(); % struct to keep all the individual neuron_a info
    % initialize the struct for this neuron
    neuron_a.valley_end =[];
    neuron_a.valley_start = [];
    neuron_a.plateau_start = [];
    neuron_a.plateau_end = [];
    neuron_a.valley_start_c = [];
    neuron_a.valley_end_c = [];

    if y(1)==0 % if the 1st frame is 0 assign it to valley start
        neuron_a.valley_start = [neuron_a.valley_start; 1];
    end

    for i=2:length(y)-1 
        % if current frame is >0 and previous is 0, valley ended in i-1
        if y(i)>0 && y(i-1)==0
            neuron_a.valley_end = [neuron_a.valley_end; i-1];
            previous=1;
        % if current frame is 0 and previous is >0,valley starts at i
        elseif y(i)==0 && y(i-1)>0
            neuron_a.valley_start = [neuron_a.valley_start; i];
            previous=0;
        end    
    end
    
    if previous==0 % if we finish with a valley
        neuron_a.valley_end=[neuron_a.valley_end; length(y)];%267
    end
       
    % drop valleys of small duration (smaller than for_concat frames)
    valley_dur = [];
    valley_dur = neuron_a.valley_end - neuron_a.valley_start;%268
    temp = find(valley_dur<=for_concat);
    neuron_a.valley_start(temp)=[];
    neuron_a.valley_end(temp)=[];
    if (~isempty(neuron_a.valley_end) || ~isempty(neuron_a.valley_start)) 
    % use valley_start_c and valley_end_c for the next steps
     neuron_a.valley_start_c=neuron_a.valley_start(1);
     neuron_a.valley_end_c=neuron_a.valley_end(1);
     counter = 1;

     % if the interval between 2 consecutive valleys is small, which is a
     % plateau:
     for i=2:length(neuron_a.valley_start)
         temp = neuron_a.valley_start(i)-neuron_a.valley_end(i-1);
            if temp<=for_concat % if less than for_concat frames then do the concat
            %with the previous valley
            % this means that the end of the valley will be updated
            neuron_a.valley_end_c(counter) = neuron_a.valley_end(i);
        else%if we do not proceed with concatenation, end the previous valley 
            %and start the next one
            neuron_a.valley_start_c=[neuron_a.valley_start_c; neuron_a.valley_start(i)];
            neuron_a.valley_end_c=[neuron_a.valley_end_c; neuron_a.valley_end(i)];
            counter=counter+1;
        end
    end
    
    flag = 1; % flag to consider (gets 0 when 1st frame belongs to valley) or plateau
    % define plateaus 
    if neuron_a.valley_start_c(1)>1 % if the 1st frame does not belong to valley....
        neuron_a.plateau_start=[neuron_a.plateau_start; 1];
    else
        neuron_a.plateau_start=[neuron_a.plateau_start; neuron_a.valley_end_c(1)+1];
        flag = 0;
    end
    if flag==1
        for i=1:length(neuron_a.valley_start_c) % define intermediate plateaus
            neuron_a.plateau_end=[neuron_a.plateau_end; neuron_a.valley_start_c(i)-1];
            neuron_a.plateau_start=[neuron_a.plateau_start; neuron_a.valley_end_c(i)+1];
        end
    else
        for i=2:length(neuron_a.valley_start_c) % define intermediate plateaus
            neuron_a.plateau_end=[neuron_a.plateau_end; neuron_a.valley_start_c(i)-1];
            neuron_a.plateau_start=[neuron_a.plateau_start; neuron_a.valley_end_c(i)+1];
        end
    end
    neuron_a.plateau_start(end)=[];%neuron_V(n).plateau_start(end)=[];
    %     if neuron_a.valley_end_c(end)==length(y) % define the end of the last plateau
%         neuron_a.plateau_end=[neuron_a.plateau_end; neuron_a.valley_start_c+1];
%     else
%         neuron_a.plateau_end=[neuron_a.plateau_end; length(y)];
%     end
    end

