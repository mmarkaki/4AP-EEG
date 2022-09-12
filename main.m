#start from raw fluoresence files, where the number of frames per segment has been stored in a vector for ease of processing:
nframes_seg=[8000,8000,8001,8001,32004,32004,32004,32004,32004];

#also the framePeriod stands for the duration of a calcium frame in seconds, it differs per mouse and per segment; these values also
#can be read from xml relevant files and stored in vectors to load:

#Relevant preprocessing includes putting the segments in their correct order for mice 4 JL65, 5 JL66, JL06VIP and JL07VIP.
#Also mice 12 and 13 need a proper selection of cells out of provided ROIs; here, we suppose that these preprocessing steps have been
#already performed hence we can simply load the corrected raw fluoresence files and proceed in the same way for all mice automatically.

#Load raw fluoresence file f, fV:

[nframes,cellnum]=size(f);
[nframesV,cellnumV]=size(fV);

#df/f from raw fluoresence files is called dff; it is estimated per segment:




# z scored df/f is called newdff:

#Definition of threshold thresh from vehicle:

#Local plateaus definition based on  z scored df/f and  threshold defined from vehicle as:
nframes=nsamples
forConcat=100;
dffV=newdff;
[nframes,cellnum]=size(newdff);
filtV = zeros(nframes,cellnum);
smtlbV = filtV;
for n = 1:cellnum
    smtlbV(:,n) = dffV(:,n);
    filtV(:,n)=~(smtlbV(:,n)<thresh(end));%
end%
neuron_a = struct(); % struct to keep all the individual neuron_a info
for n=cellnum:-1:1
    neuron_a(n).valley_end =[];
    neuron_a(n).valley_start = [];
    neuron_a(n).plateau_start = [];
    neuron_a(n).plateau_end = [];
    neuron_a(n).valley_start_c = [];
    neuron_a(n).valley_end_c =[]; 
end%"problems" are due to that neurons are below these thresholds
for n=1:cellnum%
    y = filtV(1:nframes,n).*(smtlbV(1:nframes,n));
    if isempty(find(y, 1))
        continue;
    end
    neuron_a(n) = local_plateaus_neuron(y,forConcat);
    %when valleys shorter than for_concat, they are all
    %dropped and plateaus are formed
end
%statistics on the number and duration of plateaus:
numb=zeros(1,cellnum);%maxdur=zeros(1,cellnumV);
for n=1:cellnum%length(neuron_V(n).plateau_end)
    y = filtV(1:nframes,n).*(smtlbV(1:nframes,n));
    if isempty(find(y, 1))
        continue;
    end
    neuron_a(n).duration = neuron_a(n).plateau_end-neuron_a(n).plateau_start;%maxdur(n)=max(neuron_V(n).duration);%2334 frames
    numb(n)=length(neuron_a(n).plateau_start);%how many plateus
end
