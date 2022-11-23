clear all
close all
mouse_id=7;
k=1
forConcatP=100;
forConcatV=100;
load(sprintf('E:/matfiles/mouse%d/4AP/neuronbLocalPlateaus.mat',mouse_id))
load(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_start_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k))%,'plateau_start')
load(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_end_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k))%,'plateau_end')

neuron_med=neuron;
cellnum=length(neuron_med);
clear neuron
onsets=[];
offsets=[];
Rsquares_med=[];
glob_med=[];
Amplitds_med=[];
lag_med=[];%Rmax/2 / slope
lag_glob_med=[];
lag2_med=[];%X50_med-(onset-start)
SlopeDeriv_med=[];
neuron_rec_med=[];
neuron_rec_med2=[];
plat_rec_med=[];
rsqS=0.8

for n=1:cellnum%separate the properly fitted values based on Rsquare first:
         j=find(neuron_med(n).rsqSigm>=rsqS);%140
         neuron_rec_med=[neuron_rec_med length(neuron_med(n).start(j))];
         for jj=1:length(j)
             neuron_rec_med2=[neuron_rec_med2 n];%neuron n appears length(j) times
         end
         plat_rec_med=[plat_rec_med neuron_med(n).start(j)];
         Rsquares_med=[Rsquares_med neuron_med(n).rsqSigm(j)];
         glob_med=[glob_med neuron_med(n).glob(j)];%the global plateaus
         Amplitds_med=[Amplitds_med neuron_med(n).ZmaxSG(j)];
         Rmax=neuron_med(n).ZmaxSG(j);
         temp= tan(neuron_med(n).slopeDeriv(j));%df/f z score per frames
         SlopeDeriv_med=[SlopeDeriv_med temp]; 
         lag=(Rmax/2)./temp;%slope;%
         lag_med=[lag_med lag];%Rmax/2 / slope gives one option of x50
         x50=neuron_med(n).x50True(j);
         onset=neuron_med(n).plateau_start(j)';%I need this
         offset=neuron_med(n).plateau_end(j)';%I need this too
         onsets=[onsets onset];%starts
         offsets=[offsets offset];
         start=neuron_med(n).start(j);
         lag2=x50-(onset-start);%this is more correct?
         lag2_med=[lag2_med lag2./33];%
         lag_glob_med=[lag_glob_med x50+start];%absolute time to half Rmax;
end

%FOV plots of parameters in bursts:
%Walden figs mouse 4 sigmoid, 15/10/2021:
%the number of global plateaus
nbursts=max(glob_med)%8 bursts mouse 5; 10 bursts mouse 6
%To detect cells with a reliable recruitment pattern across multiple ictal
%events, we took two different approaches. The first approach involved the
%division of each optical seizure break-in of an experiment into three equal 
%time bins (early, intermediate and late, or 1, 2 and 3, Fig. S2 C, top panel),
%yielding a single “bin” per cell per seizure event. Similarly to a z
% score, we also calculated a relative recruitment score (rr-score) for every
%single cell as follows.% With respect to the median frame within each ictal 
%onset period (Y50, the frame wherein the cumulative number of recruited cells
%first equaled or exceeded 50% of all analyzed cells), each cell was assigned
%an onset frame lag (Xn = individual cell onset frame. Frame lag = Xn – Y50).
% Each cell’s lag was then divided by the standard deviation of all lags across
%all active cells within a seizure event (?). Specifically, rr-scores were calculated 
%as Z = (Xn – Y50) / ?. From the resulting “observed recruitment matrices” (e.g.
%100x10 matrix in the case of 100 cells and 10 ictal events) of either time 
%bin categories (1, 2, or 3) or rr-scores, the standard deviation (std) was 
%obtained for every individual cell across all events (resulting in a 100x1 
%vector).
tot_neurons=[];
bigstart=[];
biglagC=[];
burstID=[];
neuronsID=[];
onsetsCor=[];
Y50s=[];
%offset4=[];
figure
hold on
for ib=1:nbursts
    jj=find(glob_med==ib);%
    tot_neurons=[tot_neurons length(jj)];
    %load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));
    if (length(jj)>3)%at least 4 neurons with local plateaus in this burst
        %in a composite burst, which are many unimodal bursts very close together,
        %a neuron might be recruited later in the composite burst, 
        %but first in a later unimodal burst
        temp=neuron_rec_med2(jj);%these are the ON neurons
        neuronsID=[neuronsID temp];
        %clag=offsets(jj)./30;%disengagement
        clag=onsets(jj)./30;%engagement
        onsetsCor=[onsetsCor onsets(jj)];
        clag(clag<0)=0;
        refere=median(clag);
        Y50s=[Y50s 30*refere*ones(1,length(jj))];%glob meds
        %stde=std(clag);%within a burst-relative disengagement
        %stde=max(abs(clag-refere));  
        clagZ=(clag-refere);%./stde;%rr-ratio of neurons in this burst
        biglagC=[biglagC clagZ];%lags values
        burstID=[burstID ib*ones(1,length(jj))];%glob meds
        %offset7=[offset7 clagZ];%4288
        if length(clagZ)>5
            [y,x]=ecdf(clagZ);
            plot(x,100*y,'Color',[1, 0, 0, 0.2])%red for recruitment
            %plot(x,100*y,'Color', [0, 1, 0, 0.2])%blue for disengagegement
        end
    end
end
ylabel('ECDF')%sprintf('
title(sprintf('Mouse %d',mouse_id))
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/ecdfLags%d.png',mouse_id,mouse_id))
         
save(sprintf('matfiles/mouse%d/4AP/mouse%d_iblagC.mat',mouse_id,mouse_id),'biglagC')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_burstid.mat',mouse_id,mouse_id),'burstID')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_neuronsid.mat',mouse_id,mouse_id),'neuronsID')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_onsets.mat',mouse_id,mouse_id),'onsetsCor')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_Y50sec.mat',mouse_id,mouse_id),'Y50s')
