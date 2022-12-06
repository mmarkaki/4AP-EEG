#start from raw fluorescence files, where the number of frames per segment has been stored in a vector for ease of processing:
nframes_seg=[8000,8000,8001,8001,32004,32004,32004,32004,32004];

#also the framePeriod stands for the duration of a calcium frame in seconds, it differs per mouse and per segment; these values
#can be read from relevant xml nd stored in vectors to load:
%Mouse 5:framePeriods=[0.033709258,0.033716781,0.033716173,0.033716891,0.033717511,0.033717443,0.033718333,0.033718122,0.033717802,0.033717893,0.033718102,0.033718176,0.033705027];
load(sprintf('E:/matfiles/mouse%d/4AP/framePeriods.mat',mouse_id))

#coordinates:
load(sprintf('E:/data/mouse%d/4AP/coordinates.mat',mouse_id))%cord=coordinates(cells4AP,:);

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

%prepare to estimate global, population plateaus/bursts based on all local plateaus, not only the ones with R squared >=0.8:
timeseriesl=zeros(size(newdff)); % here we will keep the 0-1 timeseries
for n=1:cellnum % timeseries with 0 at valleys and 1 at plateaus
    for i=1:length(neuron(n).plateau_start)
        timeseriesl(neuron(n).plateau_start(i):neuron(n).plateau_end(i),n)=1;
    end
end
%plot this:
figure('DefaultAxesFontSize',14)
plot(100*sum(timeseriesl,2)./cellnum)

%now estimate population plateaus based on timeseriesl:
num_plateaus = [];
kk=1:30%check these engagement percentages as potential cutoffs:
forConcatV=100;
forConcatP=100;%only for pop plateaus; otherwise forConcatP=211 for local plateaus
for k=1:length(kk) 
    [plateau_start, plateau_end, timeseries_gl] = global_plateaus_init(forConcatV,forConcatP,k,timeseriesl);
    num_plateaus =[num_plateaus length(plateau_start)];
    save(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_start_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k),'plateau_start')
    save(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_end_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k),'plateau_end')
end

%Choice of cutoff could be based on the maximization of the number of global plateaus in this plot: 
figure('DefaultAxesFontSize',14)
plot(kk,num_plateaus)
axis tight
ylabel('# population plateaus')
xlabel('engagement (%)')
title(sprintf('4AP mouse %d',mouse_id))
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/engagPlateausZforConcatV%d_forConcatP%d.png',mouse_id,forConcatV,forConcatP))
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/engagPlateausZforConcatV%d_forConcatP%d.fig',mouse_id,forConcatV,forConcatP))
