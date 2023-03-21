function [biglagC,burstID,neuronsID,onfsetsCor,Y50s] = ecdf_recruit_disengage(mouse_id,k,plateau_start,plateau_end,neuron,flag)
%if flag = 1 recruitment; if flag = 0 disengagement
%clear all
%close all
%mouse_id=7;
%k=1
%forConcatP=100;
%forConcatV=100;
%load(sprintf('E:/matfiles/mouse%d/4AP/neuronbLocalPlateaus.mat',mouse_id))
%load(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_start_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k))%,'plateau_start')
%load(sprintf('E:/matfiles/mouse%d/4AP/pop_plateaus_end_forConcatV%d_forConcatP%d_k%d.mat',mouse_id,forConcatV,forConcatP,k))%,'plateau_end')

cellnum=length(neuron);
onsets=[];
offsets=[];
Rsquares=[];
glob=[];
Amplitds=[];
%lag_med=[];%Rmax/2 / slope
%lag_glob_med=[];
%lag2_med=[];%X50_med-(onset-start)
SlopeDeriv=[];
neuron_rec=[];
neuron_rec2=[];
plat_rec=[];
rsqS=0.8

for n=1:cellnum%separate the properly fitted values based on Rsquare first:
         j=find(neuron(n).rsqSigm>=rsqS);%140
         neuron_rec=[neuron_rec length(neuron(n).start(j))];
         for jj=1:length(j)
             neuron_rec2=[neuron_rec2 n];%neuron n appears length(j) times
         end
         plat_rec=[plat_rec neuron(n).start(j)];
         Rsquares=[Rsquares neuron(n).rsqSigm(j)];
         glob=[glob neuron(n).glob(j)];%the global plateaus
         Amplitds=[Amplitds neuron(n).ZmaxSG(j)];
         Rmax=neuron(n).ZmaxSG(j);
         temp= tan(neuron(n).slopeDeriv(j));%df/f z score per frames
         SlopeDeriv=[SlopeDeriv temp]; 
         %lag=(Rmax/2)./temp;%slope;%
         %lag_med=[lag_med lag];%Rmax/2 / slope gives one option of x50
         %x50=neuron(n).x50True(j);
         onset=neuron(n).plateau_start(j)';%I need this
         offset=neuron(n).plateau_end(j)';%I need this too
         onsets=[onsets onset];%starts
         offsets=[offsets offset];
         %start=neuron_med(n).start(j);
         %lag2=x50-(onset-start);%this is more correct?
         %lag2_med=[lag2_med lag2./33];%
         %lag_glob_med=[lag_glob_med x50+start];%absolute time to half Rmax;
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
onfsetsCor=[];
Y50s=[];
figure
hold on
for ib=1:nbursts
    jj=find(glob==ib);%
    tot_neurons=[tot_neurons length(jj)];
    if (length(jj)>3)%at least 4 neurons with local plateaus in this burst
        %in a composite burst, which are many unimodal bursts very close together,
        %a neuron might be recruited later in the composite burst, 
        %but first in a later unimodal burst
        temp=neuron_rec2(jj);%these are the ON neurons
        neuronsID=[neuronsID temp];
        if flag==1
            clag=onsets(jj);%engagement
         else
            clag=offsets(jj);%disengagement
         end
         onfsetsCor=[onfsetsCor clag];
        refere=median(clag);
        Y50s=[Y50s refere*ones(1,length(jj))];%glob meds
        %stde=std(clag);%within a burst-relative disengagement
        %stde=max(abs(clag-refere));  
        clagZ=(clag-refere);%./stde;%rr-ratio of neurons in this burst
        biglagC=[biglagC clagZ];%lags values
        burstID=[burstID ib*ones(1,length(jj))];%glob meds
        if length(clagZ)>5
            [y,x]=ecdf(clagZ);
            if flag==1
               plot(x,100*y,'Color',[1, 0, 0, 0.2])%red for recruitment
            else
               plot(x,100*y,'Color', [0, 1, 0, 0.2])%blue for disengagegement
        end
    end
end
ylabel('ECDF')%sprintf('
title(sprintf('Mouse %d',mouse_id))
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/ecdfLags%d.png',mouse_id,mouse_id))
         
save(sprintf('matfiles/mouse%d/4AP/mouse%d_iblagC.mat',mouse_id,mouse_id),'biglagC')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_burstid.mat',mouse_id,mouse_id),'burstID')
save(sprintf('matfiles/mouse%d/4AP/mouse%d_neuronsid.mat',mouse_id,mouse_id),'neuronsID')
if flag==1
    save(sprintf('matfiles/mouse%d/4AP/mouse%d_onsets.mat',mouse_id,mouse_id),'onfsetsCor')
else
    save(sprintf('matfiles/mouse%d/4AP/mouse%d_offsets.mat',mouse_id,mouse_id),'onfsetsCor')
end
save(sprintf('matfiles/mouse%d/4AP/mouse%d_Y50sec.mat',mouse_id,mouse_id),'Y50s')
end

%Proceed to anova tests of statistical significance of spatiotemporal clustering:
save('results/mouse7/Bursts/neuronInfo.mat','neuron')
load('results/mouse7/Bursts/neuronInfo.mat')
load('../data/mouse7/4AP/coordinates.mat')
cellnum=length(neuron)
%info for onsets:
%neuron.plateaus_per_burst=[];
%Figure 2B
%ecdf init vs prop area
for ib = 1:14
    iblag=[];
    near_neurons=[];
    for i=1:cellnum
%         if (coordinates(i,1)>= 450 && coordinates(i,2)>=150 )%init area
%             continue
%         end
        ibindex=find(neuron(i).firstburst==ib);%8,ibindex=1
        if (isempty(ibindex))
             continue
         end
         near_neurons = [near_neurons i];
         iblag=[iblag neuron(i).firstlag(ibindex)];%311 - ecdf of lags for burst 8
    end
    save(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib),'iblag')
    save(sprintf('results/mouse7/Bursts/burst%d_nearn.mat',ib),'near_neurons')
    %save(sprintf('burst%d_lagsInit.mat',ib),'iblag')
    %save(sprintf('burst%d_nearnInit.mat',ib),'near_neurons')
end
figure('DefaultAxesFontSize',15)
hold on
for ib=1:14
    load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));
    refere=median(iblag);
    %stde=6*iqr(iblag);
    stde=max(abs(iblag-refere));
    iblagC=(iblag-refere)./stde;
    [y,x]=ecdf(iblagC);
    plot(x,100*y,'Color', [8 8 8]/255)
    pause
end
title('initiation area')
title('All FOV')
title('propagation area')

xlabel('sec')
ylabel('% recruited cells')
xlabel('> 50% cells recruited')
save('figures/mouse7/recruitment/wenzelb2.fig')
save('figures/mouse7/recruitment/wenzelb2.png')
save('figures/mouse5/recruitment/wenzelb2Init.fig')
save('figures/mouse5/recruitment/wenzelb2Init.png')
save('figures/mouse5/recruitment/wenzelb2Prop.fig')
save('figures/mouse5/recruitment/wenzelb2Prop.png')

%Figures 2E 2F
%Determination of recruitment durations was done by calculating the time 
%period from the first to last recruited registered cell, excluding the 5%
%most deviant cells.
figure('DefaultAxesFontSize',15)
hold on 
difburst=zeros(1,8);
for ib=1:14
    load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));%89
    refere1=quantile(iblag,0.05);%78.95
    refere2=quantile(iblag,0.95);%563
    newblag=iblag(intersect(find(iblag>=refere1), find(iblag<=refere2)));%81
    difburst(ib)=(max(newblag)-min(newblag))/30;
end
sum(difburst)/8%59.7/61.4/37.9/24.6917
std(difburst)%65.2/64.6/51.4/17.7065
mydat=[difburst;4*ones(1,14)];
plot(mydat(2,:),mydat(1,:),'bo')
figure(20);hold on;
plot(difburst,'bo')
legend('FOV','init.','prop.')
legend('mouse 4','mouse 5')
title('initiation area')
title('All FOV')
title('propagation area')
ylabel('recruitment duration (sec)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
%xlabel('bursts')

figure('DefaultAxesFontSize',15)
plot(difburst,'ro')
xlabel('Bursts')
ylabel('recruitment duration (sec)')
save('figures/mouse6/recruitment/wenzele2.fig')
save('figures/mouse6/recruitment/wenzele2.png')
save('figures/mouse4/wenzelb2Prop.fig')
save('figures/mouse4/wenzelb2Prop.png')

%Figures 2G 2H
%Left, Spatial analysis of propagation area: 
%Spatiotemporal quartile clustering (quartiles calculated as mean coordinate
%of 1–25% earliest cells, 25–50%, 50–75% and 75–100%, see Materials and Methods)
%across 8 consecutive seizures (bivariate ANOVA p < 0.001, all extrafocal 
figure('DefaultAxesFontSize',15)
hold on 
X1=[];
X2=[];
X1L=[];
X2L=[];
for ib=1:14
    load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));%367 lag
    load(sprintf('results/mouse7/Bursts/burst%d_nearn.mat',ib));%503 neuron
    if length(iblag)>3
        refere1=quantile(iblag,0.25);%earliest 25% 367
        newblag1=find(iblag<=refere1);%1
        location11=sum(coordinates(near_neurons(newblag1),1))/length(newblag1);
        location12=sum(coordinates(near_neurons(newblag1),2))/length(newblag1);
        plot(location12,location11,'b.','MarkerSize',20)
        X1=[X1 location11];
        X2=[X2 location12];
        refere2=quantile(iblag,0.5);%563
        newblag2=intersect(find(iblag>=refere1), find(iblag<=refere2));%81
        location11=sum(coordinates(near_neurons(newblag2),1))/length(newblag2);
        location12=sum(coordinates(near_neurons(newblag2),2))/length(newblag2);
        plot(location12,location11,'g.','MarkerSize',20)
        X1=[X1 location11];
        X2=[X2 location12];
        refere3=quantile(iblag,0.75);%563
        newblag3=intersect(find(iblag>=refere2), find(iblag<=refere3));%81
        location11=sum(coordinates(near_neurons(newblag3),1))/length(newblag3);
        location12=sum(coordinates(near_neurons(newblag3),2))/length(newblag3);
        plot(location12,location11,'y.','MarkerSize',20)
        X1L=[X1L location11];
        X2L=[X2L location12];
        %latest:
        newblag4=find(iblag>refere3);%81
        location11=sum(coordinates(near_neurons(newblag4),1))/length(newblag4);
        location12=sum(coordinates(near_neurons(newblag4),2))/length(newblag4);
        plot(location12,location11,'r.','MarkerSize',20)
        X1L=[X1L location11];
        X2L=[X2L location12];
    end
end
title('ALL FOV')
title('initiation area')
title('propagation area')
xlim(gca,[0 max(coordinates(:,2)+20)])
ylim(gca,[0 max(coordinates(:,1)+20)])
%  c.Label.String = 'Early vs Late Cells';
%    caxis%-0.3 ([62.8,639.1])16 111841
    xlabel('y-coordinates (ìm)')
    ylabel('x-coordinates (ìm)')
    legend('earliest','early','late','latest')
    save('figures/mouse7/recruitment/wenzelH2.fig')
save('figures/mouse7/recruitment/wenzelH2.png')
%Initiation area x-coordinates > 470, y coordinates > 200
%participates in all bursts coordinates(238,:)  473.3871  338.9355
%clustering of these according to Wenzel:
% X(:,1)=[X1, X1L]
% X(:,2)=[X2, X2L]
% opts = statset('Display','iter');
% [idx,C,sumd,d,midx,info] = kmedoids(X,2,'Distance','cityblock','Options',opts);

%exp. p < 0.05). Right, Spatial analysis of initiation site: Spatiotemporal 
%quartile clustering (each quartile coordinate = spatial mean of 25% recruited cells)
%clustering across 10 consecutive seizures (bivariate ANOVA p = 0.0145, all
%intrafocal exp. p < 0.05). Depiction of statistical significance for Figure 2:
%*p < 0.05; ***p < 0.001.
%a conserved spatial pattern of relative cell recruitment was evident in both
%compartments across seizures. … To quantify and compare these observed 
%spatiotemporal patterns, we used a 2-dimensional ANOVA (Wenzel et al., 2017),
%categorizing cells into temporal quartiles and comparing the variance of the
%distance of each cell with the quartile spatial mean and with the variance 
%of the distance to the spatial mean of all cells (Fig. 2G). The analysis
%yielded significant differences between distributions, with bivariate F-values
%for all experiments in the propagation area (F = 12.67, 11.64, 41.22; 
%all p < 0.001) and in the seizure initiation site. 

%distance from population centroid:
figure('DefaultAxesFontSize',15)
hold on 
X1R=[];
X2R=[];
X1Y=[];
X2Y=[];
X1G=[];
X2G=[];
X1B=[];
X2B=[];
for ib=1:14
    figure('DefaultAxesFontSize',15)
    hold on 
    load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));%235 lags
    load(sprintf('results/mouse7/Bursts/burst%d_nearn.mat',ib));%503 neuron
    if length(iblag)>3
        refere1=quantile(iblag,0.25);%earliest 25% 696
        newblag1=find(iblag<=refere1);%59
        location11=sum(coordinates(near_neurons(newblag1),1))/length(newblag1);
        location12=sum(coordinates(near_neurons(newblag1),2))/length(newblag1);
        for jj=1:length(newblag1)
            plot(coordinates(near_neurons(newblag1(jj)),2),coordinates(near_neurons(newblag1(jj)),1),'b.','MarkerSize',8)
        end
        plot(location12,location11,'b.','MarkerSize',20)
        X1B=[X1B location11];
        X2B=[X2B location12];
        refere2=quantile(iblag,0.5);%959
        newblag2=intersect(find(iblag>=refere1), find(iblag<=refere2));%59
        location11=sum(coordinates(near_neurons(newblag2),1))/length(newblag2);
        location12=sum(coordinates(near_neurons(newblag2),2))/length(newblag2);
        plot(location12,location11,'g.','MarkerSize',20)
        X1G=[X1G location11];
        X2G=[X2G location12];
        for jj=1:length(newblag2)
            plot(coordinates(near_neurons(newblag2(jj)),2),coordinates(near_neurons(newblag2(jj)),1),'g.','MarkerSize',8)
        end
        refere3=quantile(iblag,0.75);%1164
        newblag3=intersect(find(iblag>=refere2), find(iblag<=refere3));%61
        location11=sum(coordinates(near_neurons(newblag3),1))/length(newblag3);
        location12=sum(coordinates(near_neurons(newblag3),2))/length(newblag3);
        plot(location12,location11,'y.','MarkerSize',20)
        for jj=1:length(newblag3)
            plot(coordinates(near_neurons(newblag3(jj)),2),coordinates(near_neurons(newblag3(jj)),1),'y.','MarkerSize',8)
        end
        X1Y=[X1Y location11];
        X2Y=[X2Y location12];
        %latest:
        newblag4=find(iblag>refere3);%57
        location11=sum(coordinates(near_neurons(newblag4),1))/length(newblag4);
        location12=sum(coordinates(near_neurons(newblag4),2))/length(newblag4);
        plot(location12,location11,'r.','MarkerSize',20)
        X1R=[X1R location11];
        X2R=[X2R location12];
        for jj=1:length(newblag4)
            plot(coordinates(near_neurons(newblag4(jj)),2),coordinates(near_neurons(newblag4(jj)),1),'r.','MarkerSize',8)
        end
    end
    title(sprintf('Burst %d',ib));%     set(gca,'ydir','reverse')
    xlim(gca,[0 max(coordinates(:,2)+20)])
    ylim(gca,[0 max(coordinates(:,1)+20)])
    xlabel('y-coordinates (ìm)')
    ylabel('x-coordinates (ìm)')
    pause
    saveas(gcf,sprintf('figures/mouse7/recruitment/neuronsBurstQuantLags%d.png',ib))
    saveas(gcf,sprintf('figures/mouse7/recruitment/neuronsBurstQuantLags%d.fig',ib))
end
title('ALL FOV')
% title('initiation area')
% title('propagation area')
set(gca,'ydir','reverse')
xlim(gca,[0 max(coordinates(:,2)+20)])
ylim(gca,[0 max(coordinates(:,1)+20)])
%  c.Label.String = 'Early vs Late Cells';
%    caxis%-0.3 ([62.8,639.1])16 111841
xlabel('y-coordinates (ìm)')
ylabel('x-coordinates (ìm)')
legend('earliest','early','late','latest')
save('figures/mouse7/recruitment/wenzel2019H2.fig')
save('figures/mouse7/recruitment/wenzel2019H2.png')

%2-dimensional ANOVA (Wenzel et al., 2017),
%categorizing cells into temporal quartiles and comparing the variance of the
%distance of each cell with the quartile spatial mean (X1B,X1G,X1Y,X1R) and with the variance 
%of the distance to the spatial mean of all cells (mean1cellsBurst) (Fig. 2G). The analysis
%yielded significant differences between distributions, with bivariate F-values
%for all experiments in the propagation area (F = 12.67, 11.64, 41.22; 
%all p < 0.001) and in the seizure initiation site. 
%spatial mean of all cells:
% mean1cells=sum(coordinates(:,1))/cellnum;%189.59
% mean2cells=sum(coordinates(:,2))/cellnum;%304.75
mean1cellsBurst=0.25*(X1B+X1G+X1Y+X1R);%the mean x coordinate in every burst
mean2cellsBurst=0.25*(X2B+X2G+X2Y+X2R);%the mean y coordinate in every burst

bigDistB=[];
 bigDistBU=[];
  bigGroupB=[];
  bigDistG=[];
 bigDistGU=[];
  bigGroupG=[];
  bigDistY=[];
 bigDistYU=[];
  bigGroupY=[];
  bigDistR=[];
 bigDistRU=[];
  bigGroupR=[];
for ib=1:14%
    %ib=16;
    load(sprintf('results/mouse7/Bursts/burst%d_lags.mat',ib));
    %235/15/23/56/45/4/23/6/18/5/12/27/14/5/103 lags
    load(sprintf('results/mouse7/Bursts/burst%d_nearn.mat',ib));%
    if length(iblag)>3
        refere1=quantile(iblag,0.25);%earliest 25% 
        newblag1=find(iblag<=refere1);%59/4/6/14/11/3/6/2/5/1/3/7/4/1/26
        location11=sum(coordinates(near_neurons(newblag1),1))/length(newblag1);
        location12=sum(coordinates(near_neurons(newblag1),2))/length(newblag1);
        %distances from clusters centroids:
        distCellB=zeros(1,length(newblag1));
        distCellBuni=zeros(1,length(newblag1));
        %ib=15
        for i=1:length(newblag1)
            distCellB(i)=sqrt((coordinates(near_neurons(newblag1(i)),1)-location11).^2+(coordinates(near_neurons(newblag1(i)),2)-location12).^2);
            distCellBuni(i)=sqrt((coordinates(near_neurons(newblag1(i)),1)-mean1cellsBurst(ib)).^2+(coordinates(near_neurons(newblag1(i)),2)-mean2cellsBurst(ib)).^2);
        end
        bigDistB=[bigDistB distCellB];
        bigDistBU=[bigDistBU distCellBuni];
        bigGroupB=[bigGroupB length(newblag1)];
        meandistCellB(ib)=mean(distCellB);%
        meandistCellBuni(ib)=mean(distCellBuni);%
        stddistCellB(ib)=std(distCellB);%
        stddistCellBuni(ib)=std(distCellBuni);%
        refere2=quantile(iblag,0.5);%959
        newblag2=intersect(find(iblag>=refere1), find(iblag<=refere2));%59
        location11=sum(coordinates(near_neurons(newblag2),1))/length(newblag2);
        location12=sum(coordinates(near_neurons(newblag2),2))/length(newblag2);
        for i=1:length(newblag2)
            distCellG(i)=sqrt((coordinates(near_neurons(newblag2(i)),1)-location11).^2+(coordinates(near_neurons(newblag2(i)),2)-location12).^2);
            distCellGuni(i)=sqrt((coordinates(near_neurons(newblag2(i)),1)-mean1cellsBurst(ib)).^2+(coordinates(near_neurons(newblag2(i)),2)-mean2cellsBurst(ib)).^2);
        end
        bigDistG=[bigDistG distCellG];
        bigDistGU=[bigDistGU distCellGuni];
        bigGroupG=[bigGroupG length(newblag2)];
        
        meandistCellG(ib)=mean(distCellG);%%
        meandistCellGuni(ib)=mean(distCellGuni);%
        stddistCellG(ib)=std(distCellG);%
        stddistCellGuni(ib)=std(distCellGuni);%
        refere3=quantile(iblag,0.75);%1164
        newblag3=intersect(find(iblag>=refere2), find(iblag<=refere3));%61
        location11=sum(coordinates(near_neurons(newblag3),1))/length(newblag3);
        location12=sum(coordinates(near_neurons(newblag3),2))/length(newblag3);
        for i=1:length(newblag3)
            distCellY(i)=sqrt((coordinates(near_neurons(newblag3(i)),1)-location11).^2+(coordinates(near_neurons(newblag3(i)),2)-location12).^2);
            distCellYuni(i)=sqrt((coordinates(near_neurons(newblag3(i)),1)-mean1cellsBurst(ib)).^2+(coordinates(near_neurons(newblag3(i)),2)-mean2cellsBurst(ib)).^2);
        end
        bigDistY=[bigDistY distCellY];
        bigDistYU=[bigDistYU distCellYuni];
        bigGroupY=[bigGroupY length(newblag3)];
         meandistCellY(ib)=   mean(distCellY);%
        meandistCellYuni(ib)=mean(distCellYuni);%
        stddistCellY(ib)=std(distCellY);%
        stddistCellYuni(ib)=std(distCellYuni);%
        %latest:
        newblag4=find(iblag>refere3);%57
        location11=sum(coordinates(near_neurons(newblag4),1))/length(newblag4);
        location12=sum(coordinates(near_neurons(newblag4),2))/length(newblag4);
        for i=1:length(newblag4)
            distCellR(i)=sqrt((coordinates(near_neurons(newblag4(i)),1)-location11).^2+(coordinates(near_neurons(newblag4(i)),2)-location12).^2);
            distCellRuni(i)=sqrt((coordinates(near_neurons(newblag4(i)),1)-mean1cellsBurst(ib)).^2+(coordinates(near_neurons(newblag4(i)),2)-mean2cellsBurst(ib)).^2);
        end
        bigDistR=[bigDistR distCellR];
        bigDistRU=[bigDistRU distCellRuni];
        bigGroupR=[bigGroupR length(newblag4)];
        meandistCellR(ib)=mean(distCellR);%
        meandistCellRuni(ib)=mean(distCellRuni);%
        stddistCellR(ib)=std(distCellR);%
        stddistCellRuni(ib)=std(distCellRuni);%
    end
end

%categorizing cells into temporal quartiles and comparing the variance of 
%the distance of each cell with the quartile spatial mean and with the 
%variance of the distance to the spatial mean of all cells 
%anova%y1=[stddistCellB,stddistCellG,stddistCellY,stddistCellR]
y1=[bigDistB,bigDistG,bigDistY,bigDistR]%2805
y2=[bigDistBU,bigDistGU,bigDistYU,bigDistRU]%2805%y2=[stddistCellBuni,stddistCellGuni,stddistCellYuni,stddistCellRuni]
y=[y1;y2]
%y=[y2;y1];
%group=
group=[ones(1,105),2*ones(1,231),3*ones(1,224),4*ones(1,249)];
%group=[ones(1,15),2*ones(1,15),3*ones(1,15),4*ones(1,15)];
[P,ANOVATAB,STATS] = anova1(y,group)
[c,m,h,nms] = multcompare(STATS,'display','off');
 [nms num2cell(m)]   
%number of members:
n=[235,15,23,56,45,4,23,6,18,5,12,27,14,5,103]
%anova in each burst separately:
y1=[bigDistB]%2805
y2=[bigDistBU]%2805%y2=[stddistCellBuni,stddistCellGuni,stddistCellYuni,stddistCellRuni]
y=[y1;y2]
%y=[y2;y1];
%group=
group=[ones(1,337)];
%group=[ones(1,15),2*ones(1,15),3*ones(1,15),4*ones(1,15)];
[P,ANOVATAB,STATS] = anova1(y,group)
[c,m,h,nms] = multcompare(STATS,'display','off');
 [nms num2cell(m)]   
 


