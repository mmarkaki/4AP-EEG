%framePeriods for calcium provided in xml files; stored in mat files;
load(sprintf('E:/matfiles/mouse%d/4AP/framePeriods.mat',mouse_id))
seg_framesCa=8000;%seg_framesCa=8001;seg_length=32004;%seg_length=8000;seg_length=8001;
CaSeg=seg_framesCa*framePeriods%about 269.6 sec
load(sprintf('E:/data/mouse%d/4AP/Tseries_VoltageRecording_001.mat',mouse_id))
eeg4AP=eegData;
clear eegData
Fs = eeg4AP.Fs;
seg_frames=length(eeg4AP.CleanEEG_cV(:,1))/seg_num
EEGsegSec=seg_frames/Fs
minfreq=1;
maxfreq=100;
Fc = 60;%noise to be filtered out from EEG
Wc = Fc/(Fs/2);
BW = Wc/10;%35;
[b,a] = iirnotch(Wc(1),BW(1));%the filter characteristics
%we initially consider a whole segment
%then, we could choose windows within this segment, depending on the application, eg, 50 sec sliding non-overlapping windows:
fiftysec=50*Fs;%EEG frames
fiftysecCa=floor(50./framePeriods)%Ca frames
nwind=floor(seg_frames/fiftysec)%number of 50 sec sliding non-overlapping windows
for iseg=1:seg_num
    tstart0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
    tend0=iseg*seg_frames;
    eegTimeSegment = eeg4AP.Time_s(tstart0:tend0);
    yinitial=eeg4AP.CleanEEG_cV(tstart0:tend0, 1);%whole segment
    for j=1:nwind
        tstart=1 + (j-1)*fiftysec;
        t_end=j*fiftysec;
        x1=yinitial(tstart:t_end,1);
        time1=eegTimeSegment(tstart:t_end);
        tstartCa=1+(j-1)*fiftysecCa(iseg)+(iseg-1)*seg_length;%corresponding calcium start
        t_endCa=j*fiftysecCa(iseg)+(iseg-1)*seg_length;%end
        ttt=tstartCa:t_endCa;%calcium frames to consider  
       y1 = filter(b,a,x1);
        %just once compute filterbank fb for speed;Continuous Wavelet Transform of EEG segment using Morlet wavelet
        if iseg==1 & j==1
            fb = cwtfilterbank('SignalLength',length(y1),'SamplingFrequency',Fs,'FrequencyLimits',[minfreq maxfreq], 'Wavelet',"amor");
        end
        [cfs, frqs] = wt(fb,y1);%
        nfreqs=length(frqs);
        interp_frqs = log10(frqs);%no interpolation, it takes too much time
         interp_cfs =log10(abs(cfs)); %same
        figure('DefaultAxesFontSize',14)
        hold on
        subplot(3,1,1)
        plot(ttt,zcorrectC(ttt));%calcium z score
        axis tight;
        title(sprintf('Ca++ signal mouse 4AP JL%d segment %d',mouse_id,iseg'))
        ylabel('mean z score')
        xlabel('Ca++ frames')
        subplot(3,1,2)
        plot(time1,y1)
        axis tight;
        title(sprintf('EEG signal'))
        ylim([-0.15 0.15])
        ylabel('cV')
        subplot(3,1,3)
        imagesc(time1,interp_frqs, interp_cfs);
        set(gca,'YDir','normal') %cmocean('thermal') 
        set(gca,'XTickMode','auto')
        set(gca,'YTickMode','auto')
        set(gca,'TickDir','out')%%c
        caxis([-4   -1.6])%to be fine tuned for clarity of figure
        shading interp; % interpolate shading to make color transitions visually smooth
        axis tight;
        title(sprintf('CWT mouse 4AP JL%d segment %d',mouse_id,iseg))
        ylabel('Frequency (Hz)')
        yticks([0 0.5 1 1.5 2])
        yticklabels({'1','3.2','10','31.6','100'})
        xlabel('Time (sec)')    %colorbar
        saveas(gcf,sprintf('figures/mouse%d/4AP/EEG/triples50sec/CWT_seg%d_wind%d.png',mouse_id,iseg,j));
        close
    end
end 

%PSD: Estimate power spectral density via the multitaper method using Slepian or sine tapers.
%PSD everywhere, the whole EEG:17 sec in every segment of 18 minutes means 
% 64 sliding windows of 20 sec with 3 sec overlap:
half_sec=0.5*Fs
nwind = 64
n=30000;%3 sec sliding window; overlap by 2.5 sec or 25000 frames
p=25000;
bin=n-p%half a sec in EEG frames
binCa=(bin/Fs)./framePeriods%half a sec in calcium frames
half_sec=n-p%5000 EEG frames per half sec
nseg=seg_num
for iseg=1:nseg 
    tstartEEG0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
    tendEEG0=iseg*seg_frames%
    tstartCa0=1+(iseg-1)*seg_framesCa;%EEG counter  time started from zero
    tendCa0=iseg*seg_framesCa%yinitial1=eeg
    %buffers now:    for jw=1:nwind
    onsetEEG=eeg4AP.CleanEEG_cV(tstartEEG0+half_sec:tendEEG0-half_sec,1);
    [buffer1,z]=buffer(onsetEEG,n,p);%30000 *2193  windows
    buffer2=buffer1(:,6:end);%28
    nPSDs=size(buffer2,2);%just plot this, the rest is empty
    matrixOnset=zeros(nfreqs,nPSDs);% 
    confid=zeros(nfreqs,2,nPSDs);
    for metr=1:nPSDs
        signs=buffer2(:,metr);
        y1=filter(b,a,signs);
            %computes the multitaper PSD estimate at the frequencies specified in frqs which have been estimated as following:
%             fb = cwtfilterbank('SignalLength',length(y1),'SamplingFrequency',Fs,'FrequencyLimits',[minfreq maxfreq], 'Wavelet',"amor")          
%             [cfs, frqs] = wt(fb,y1);%67 frequencies
        [pxx,frqs,Pxxc]=pmtm(y1,3.5,frqs,Fs,'ConfidenceLevel',0.95);%,'freqrange');%'unity','centered');
        matrixOnset(:,metr)=pxx;%mean(interp_cfs2,2);
        confid(:,:,metr)=Pxxc;
     end
     save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/matrixOnsetSeg%d_multitaper.mat',mouse_id,iseg),'matrixOnset');
     save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/confidOnsetSeg%d_multitaper.mat',mouse_id,iseg),'confid');
     %random from the ones already estimated
end

xxx=[ceil(plateau_start/seg_framesCa)     ceil(plateau_end/seg_framesCa)]
%mouse_id=6;
%bin=5000;
%PSD plus random estimation onsets:
yinitial=eeg4AP.CleanEEG_cV(:,1);
for iseg=1:nseg    
    tstartEEG0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
    tendEEG0=iseg*seg_frames%
    tstartCa0=1+(iseg-1)*seg_framesCa;%EEG counter  time started from zero
    tendCa0=iseg*seg_framesCa;
    lips=find(xxx(:,1)==iseg);
    for ipp=1:length(lips)
        ip=lips(ipp);%ip2=lips2(ipp);
        temp=0;
        if iseg>1
            temp=EEGsegSec*(iseg-1)-sum(CaSeg(1:iseg-1));
        end
        tstart=temp*Fs+plateau_start(ip)*framePeriods(iseg)*Fs;%337);   
        tstartEEG=round(tstart)-20*bin+1;%bin is half a second
        tstartCa=round(plateau_start(ip))-round(20*binCa(iseg))+1;
        if tstartCa<0 %then also tstartEEG<0
            tstartCa=tstartCa0;%segment start
            tstartEEG=tstartEEG0;
        end    
        t_endEEG=round(tstart+20*bin);
        t_endCa=round(plateau_start(ip)+20*binCa(iseg));
        %13th and 23rd burst end at end of segments
        if t_endCa > tendCa0
            t_endCa=tendCa0;%segment end
            t_endEEG=tendEEG0;
        end
        %here newer code loads results:
        onsetEEG=eeg4AP.CleanEEG_cV(tstartEEG:t_endEEG,1);
        %load('E:/matfiles/mouse13/4AP/zscoreMeanSJL13Cells.mat')
        %onsetCa=zcorrectC(tstartCa:t_endCa);%498
        [buffer1,z]=buffer(onsetEEG,n,p);%30000x33, first 5 contain zeros; omit them;
        buffer2=buffer1(:,6:end);%28
        nPSDs=size(buffer2,2);%48/35- just plot this, the rest is empty
        matrixOnset=zeros(nfreqs,nPSDs);% matrixOnset2=zeros(90,nPSDs);
        confid=zeros(nfreqs,2,nPSDs);
        for metr=1:nPSDs
            signs=buffer2(:,metr);
            y1=filter(b,a,signs);
            %computes the multitaper PSD estimate at the frequencies specified in f
%             fb = cwtfilterbank('SignalLength',length(y1),'SamplingFrequency',Fs,'FrequencyLimits',[minfreq maxfreq], 'Wavelet',"amor")          
%             [cfs, frqs] = wt(fb,y1);%67 frequencies
            [pxx,frqs,Pxxc]=pmtm(y1,3.5,frqs,Fs,'ConfidenceLevel',0.95);%,'freqrange');%'unity','centered');
            matrixOnset(:,metr)=pxx;%mean(interp_cfs2,2);
            confid(:,:,metr)=Pxxc;
        end
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/matrixOnsetSeg%dBurst%d_multitaper.mat',mouse_id,iseg,ip),'matrixOnset');
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/confidOnsetSeg%dBurst%d_multitaper.mat',mouse_id,iseg,ip),'confid');
        %before b
%       %circular shift for null distribution:
        r=randi([1 length(yinitial)],30,1);
        RandmatrixOnset=zeros(length(frqs),nPSDs,30);
        for mmetr=1:30
            K=r(mmetr);
            yfinal=circshift(yinitial,K);
            x1=yfinal(tstartEEG:t_endEEG);%
            z=[];
            [buffer1,z]=buffer(x1,n,p);
            %omit the first 5 
            buffer2=buffer1(:,6:end);
            nPSDs=size(buffer2,2);%28
            for metr=1:nPSDs
                signs=buffer2(:,metr);
                %signsV=buffer2V(:,metr);
                y1=filter(b,a,signs);%y1 = filter(b,a,x);       
                [pxx,frqs]=pmtm(y1,3.5,frqs,Fs);%,'freqrange');%'unity','centered');
                RandmatrixOnset(:,metr,mmetr)=pxx;
             end
        end
        meanRandmatrixOnset=mean(RandmatrixOnset,3);
        stdRandmatrixOnset=std(RandmatrixOnset,0,3);
        semRandmatrixOnset=std(RandmatrixOnset,0,3)/sqrt(length(RandmatrixOnset));
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/meanRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip),'meanRandmatrixOnset');
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/stdRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip),'stdRandmatrixOnset');
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/semRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip),'semRandmatrixOnset');
        %I could just use the same with first 4-5 plateau without estimating them for every
        %plateau; just pick them at random!
%         load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/meanRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip))
%         load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/stdRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip))
        z_score=(matrixOnset-meanRandmatrixOnset)./stdRandmatrixOnset;
        save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/zscore_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip),'z_score');
        %gather z scores and keep the distributions from distance from end:
        %figures:
        figure('DefaultAxesFontSize',14)
        hold on
        imagesc([-19:20],frqs,z_score)%of onset/offset
        ylabel('Frequency (Hz)')
%         yticks([0 1 2])
%         yticklabels({'1','10','100'})
        set(gca,'YDir','normal') 
        colorbar%'PSD (dB/Hz)'
        shading interp; % interpolate shading to make color transitions visually smooth
        title(sprintf('burst %d onset PSD',ip))
        axis tight
        saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Onsets/burst%d_PSDonset.png',mouse_id,ip));
        close
    end
end
%repeat for offsets:
 nseg=seg_num
    for iseg=1:nseg 
        tstartEEG0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
        tendEEG0=iseg*seg_frames%
        tstartCa0=1+(iseg-1)*seg_framesCa;%EEG counter  time started from zero
        tendCa0=iseg*seg_framesCa%yinitial1=eegData.CleanEEG_cV(tstart0:tend0, 1);%whole segment   
        lips=find(xxx(:,2)==iseg)% %around onset and around offset:
        for ipp=1:length(lips)%1:length(plateau_start)%segments 8,10 without bursts
            ip=lips(ipp);%ip2=lips2(ipp);
            temp=0;
            if iseg>1
                temp=EEGsegSec*(iseg-1)-sum(CaSeg(1:iseg-1))
            end
            tstart=temp*Fs+plateau_end(ip)*framePeriods(iseg)*Fs;%337);%=round(((iseg-1)*(-CaSeg+EEGseg)+plateau_end(ip))*337);    
            %=round(((iseg-1)*(-CaSeg+EEGseg)+plateau_end(ip))*337);    
            tstartEEG=round(tstart)-20*bin+1;%bin is half a second
            tstartCa=plateau_end(ip)-round(20*binCa(iseg))+1;
            if tstartCa<0 %then also tstartEEG<0
                tstartCa=tstartCa0;%segment start
                tstartEEG=tstartEEG0;
            end    
            t_endEEG=round(tstart+20*bin);
            t_endCa=plateau_end(ip)+round(20*binCa);
            %13th and 23rd burst end at end of segments
            if t_endCa > tendCa0
                t_endCa=tendCa0;%segment end
                t_endEEG=tendEEG0;
            end
            %here newer code loads results:
            offsetEEG=eeg4AP.CleanEEG_cV(tstartEEG:t_endEEG);%167737/337=497.7359
            offsetCa=zcorrectC(tstartCa:t_endCa);%498
            [buffer1,z]=buffer(offsetEEG,n,p);%30000x33, first 5 contain zeros; omit them;
            buffer2=buffer1(:,6:end);%28
            nPSDs=size(buffer2,2);%48/35- just plot this, the rest is empty
            matrixOffset=zeros(nfreqs,nPSDs);%
            confid=zeros(nfreqs,2,nPSDs);
            for metr=1:nPSDs
                signs=buffer2(:,metr);
                y1=filter(b,a,signs);
                [pxx,frqs,Pxxc]=pmtm(y1,3.5,frqs,Fs,'ConfidenceLevel',0.95);%,'freqrange');%'unity','centered');
                matrixOffset(:,metr)=pxx;%mean(interp_cfs2,2);
                confid(:,:,metr)=Pxxc;
            end
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/matrixOffsetSeg%dBurst%d_multitaper.mat',mouse_id,iseg,ip),'matrixOffset');
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/confidOffsetSeg%dBurst%d_multitaper.mat',mouse_id,iseg,ip),'confid');
            %before b
            %circular shift for null distribution was done for irrelevant pieces but that's fine:
            r=randi([1 length(yinitial)],30,1);
            RandmatrixOffset=zeros(length(frqs),nPSDs,30);
            for mmetr=1:30
                K=r(mmetr);
                yfinal=circshift(yinitial,K);
                x1=yfinal(tstartEEG:t_endEEG);%
                z=[];
                [buffer1,z]=buffer(x1,n,p);
                %omit the first 5 
                buffer2=buffer1(:,6:end);
                nPSDs=size(buffer2,2);%28
                for metr=1:nPSDs
                    signs=buffer2(:,metr);
                    y1=filter(b,a,signs);%y1 = filter(b,a,x);         
                    [pxx,frqs]=pmtm(y1,3.5,frqs,Fs);%,'freqrange');%'unity','centered');
                    RandmatrixOffset(:,metr,mmetr)=pxx;
                end
            end
            meanRandmatrixOffset=mean(RandmatrixOffset,3);
            stdRandmatrixOffset=std(RandmatrixOffset,0,3);
            semRandmatrixOffset=std(RandmatrixOffset,0,3)/sqrt(length(RandmatrixOffset));
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/meanRandom_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip),'meanRandmatrixOffset');
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/stdRandom_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip),'stdRandmatrixOffset');
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/semRandom_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip),'semRandmatrixOffset');
            %I could just use the same with first 4-5 plateau without estimating them for every
            %plateau; just pick them at random!
    %         load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/meanRandom_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip));%,'meanRandmatrixOffset');
    %         load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/stdRandom_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip));%,'stdRandmatrixOffset');
            z_score=(matrixOffset-meanRandmatrixOffset)./stdRandmatrixOffset;
            save(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Offsets/zscore_%dSeg_matrixOffset_burst%d.mat',mouse_id,iseg,ip),'z_score');
            %figures:
            figure('DefaultAxesFontSize',14)
            hold on
            imagesc([-19:20],frqs,z_score)%of onset/offset
            ylabel('Frequency (Hz)')
    %         yticks([0 1 2])
    %         yticklabels({'1','10','100'})
            set(gca,'YDir','normal') 
            colorbar%'PSD (dB/Hz)'
            shading interp; % interpolate shading to make color transitions visually smooth
            title(sprintf('burst %d offset PSD',ip))
            axis tight
            saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Offsets/burst%d_PSDoffset.png',mouse_id,ip));
            close
        end
    end
    %onsets figures loading PSD results:
    %avrgConfid=nan(nfreqs,2,35,100);
    interp_frqs = frqs
delta=find(10.^interp_frqs>=1 & 10.^interp_frqs < 4)%20 frequencies
theta=find(10.^interp_frqs>=4 & 10.^interp_frqs <8)%10 frequencies
alpha=find(10.^interp_frqs>=8 & 10.^interp_frqs <12)%6 frequencies
%beta=10.^interp_frqs(10.^interp_frqs>=12 & 10.^interp_frqs <30)%13 frequencies
beta=find(10.^interp_frqs>=15 & 10.^interp_frqs <30)%10 frequencies-not between 12 to 15 Hz?
gamma1=find(10.^interp_frqs>=30 & 10.^interp_frqs <60)%10 frequencies
gamma1_5=find(10.^interp_frqs>60 & 10.^interp_frqs <=90)%6 frequencies
gamma2=find(10.^interp_frqs>=60 & 10.^interp_frqs <100)%17 frequencies/200
%from 10^1.5786=37 Hz up to 10^1.9097=81.23 Hz, filt
for iseg=5:nseg 
%     tstartEEG0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
%     tendEEG0=iseg*seg_frames%
%     tstartCa0=1+(iseg-1)*seg_framesCa;%EEG counter  time started from zero
%     tendCa0=iseg*seg_framesCa%yinitial1=eegData.CleanEEG_cV(tstart0:tend0, 1);%whole segment   
    lips=find(xxx(:,1)==iseg)%1:6/7,8 start at first segment/2 at second
    for ipp=1:length(lips)%1:length(plateau_start)%segments 8,10 without bursts
        ip=lips(ipp);%ip2=lips2(ipp);
%         tstart=(iseg-1)*(-CaSeg(iseg)+EEGsegSec)*Fs+plateau_start(ip)*framePeriods(iseg)*Fs;%337);%=round(((iseg-1)*(-CaSeg+EEGseg)+plateau_end(ip))*337);       
%         tstartEEG=tstart-20*bin+1;%bin is half a second
%         tstartCa=plateau_start(ip)-round(20*binCa(iseg))+1;
%         if tstartCa<0 %then also tstartEEG<0
%             tstartCa=tstartCa0;%segment start
%             tstartEEG=tstartEEG0;
%         end    
%         t_endEEG=tstart+20*bin;
%         t_endCa=plateau_start(ip)+round(20*binCa);
%         %13th and 23rd burst end at end of segments
%         if t_endCa > tendCa0
%             t_endCa=tendCa0;%segment end
%             t_endEEG=tendEEG0;
%         end
        load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/meanRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip))
        %;,'meanRandVmatrixOnset');   
        load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/stdRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip));
        %,'stdRandVmatrixOnset');
        load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/semRandom_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip));
        %,'semRandVmatrixOnset');
        load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/zscore_%dSeg_matrixOnset_burst%d.mat',mouse_id,iseg,ip));
        load(sprintf('E:/matfiles/mouse%d/4AP/Buffers/OnsetsOffsets/Onsets/matrixOnsetSeg%dBurst%d_multitaper.mat',mouse_id,iseg,ip));%90x35
        avrgOnset(:,:,ip)=matrixOnset;
        avrgMean(:,:,ip)=meanRandmatrixOnset;
        avrgSEM(:,:,ip)=semRandmatrixOnset;
        avrgSTD(:,:,ip)=stdRandmatrixOnset;
        avrgZscore(:,:,ip)=z_score;
    end
end
%save and plot the new averages:

%mean over bursts:
Seiz=[28,35,36,71,72];
IIE=[1:27,29:34,37:70,73:100];
averageOnset=mean(avrgOnset,3,'omitnan');
averageOnsetSeizure=mean(avrgOnset(:,:,Seiz),3,'omitnan');%bursts 26:34
averageOnsetIIE=mean(avrgOnset(:,:,IIE),3,'omitnan');%bursts 26:34
averageZscoreOnsetSeizure=mean(avrgZscore(:,:,Seiz),3,'omitnan');%bursts 26:34
averageZscoreOnsetIIE=mean(avrgZscore(:,:,IIE),3,'omitnan');%bursts 26:34
%I am interesetd to plot z scores:
averageSEMOnsetSeizure=mean(avrgSEM(:,:,Seiz),3,'omitnan');
averageSEMOnsetIIE=mean(avrgSEM(:,:,IIE),3,'omitnan');
%median over bursts:
averageOnset=median(avrgOnset,3,'omitnan');
averageOnsetSeizure=median(avrgOnset(:,:,Seiz),3,'omitnan');%bursts 26:34
averageOnsetIIE=median(avrgOnset(:,:,IIE),3,'omitnan');%bursts 26:34
averageZscoreOnsetSeizure=median(avrgZscore(:,:,Seiz),3,'omitnan');%bursts 26:34
averageZscoreOnsetIIE=median(avrgZscore(:,:,IIE),3,'omitnan');%bursts 26:34
averageSEMOnsetSeizure=median(avrgSEM(:,:,Seiz),3,'omitnan');
averageSEMOnsetIIE=median(avrgSEM(:,:,IIE),3,'omitnan');
% averageEngageOnset=median(avrgEngage,2,'omitnan');%40 bins
% averageEngageOnsetSeizure=median(avrgEngage(:,Seiz),2,'omitnan');%40 bins
% averageEngageOnsetIIE=median(avrgEngage(:,[1:25,33,35,36]),2,'omitnan');%40 bins
% iqrEngageOnset=iqr(avrgEngage,2);%40 bins
% iqrEngageOnsetSeizure=iqr(avrgEngage(:,[26:32,34]),2);%40 bins
% iqrEngageOnsetIIE=iqr(avrgEngage(:,[1:25,33,35,36]),2);%40 bins

%averages over frequency bands:G1
averageOnsetG1=mean(averageOnset(gamma1,:),1,'omitnan');
averageOnsetSeizureG1=mean(averageOnsetSeizure(gamma1,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG1=mean(averageOnsetIIE(gamma1,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG1=mean(averageZscoreOnsetSeizure(gamma1,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG1=mean(averageZscoreOnsetIIE(gamma1,:),1,'omitnan');%bursts 26:34
%G1_5
averageOnsetG1_5=mean(averageOnset(gamma1_5,:),1,'omitnan');
averageOnsetSeizureG1_5=mean(averageOnsetSeizure(gamma1_5,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG1_5=mean(averageOnsetIIE(gamma1_5,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG1_5=mean(averageZscoreOnsetSeizure(gamma1_5,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG1_5=mean(averageZscoreOnsetIIE(gamma1_5,:),1,'omitnan');%bursts 26:34
%G2
averageOnsetG2=mean(averageOnset(gamma2,:),1,'omitnan');
averageOnsetSeizureG2=mean(averageOnsetSeizure(gamma2,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG2=mean(averageOnsetIIE(gamma2,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG2=mean(averageZscoreOnsetSeizure(gamma2,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG2=mean(averageZscoreOnsetIIE(gamma2,:),1,'omitnan');%bursts 26:34
%Alpha:
averageOnsetA=mean(averageOnset(alpha,:),1,'omitnan');
averageOnsetSeizureA=mean(averageOnsetSeizure(alpha,:),1,'omitnan');%bursts 26:34
averageOnsetIIEA=mean(averageOnsetIIE(alpha,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureA=mean(averageZscoreOnsetSeizure(alpha,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEA=mean(averageZscoreOnsetIIE(alpha,:),1,'omitnan');%bursts 26:34
%Beta:
averageOnsetB=mean(averageOnset(beta,:),1,'omitnan');
averageOnsetSeizureB=mean(averageOnsetSeizure(beta,:),1,'omitnan');%bursts 26:34
averageOnsetIIEB=mean(averageOnsetIIE(beta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureB=mean(averageZscoreOnsetSeizure(beta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEB=mean(averageZscoreOnsetIIE(beta,:),1,'omitnan');%bursts 26:34
%Delta:
averageOnsetD=mean(averageOnset(delta,:),1,'omitnan');
averageOnsetSeizureD=mean(averageOnsetSeizure(delta,:),1,'omitnan');%bursts 26:34
averageOnsetIIED=mean(averageOnsetIIE(delta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureD=mean(averageZscoreOnsetSeizure(delta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIED=mean(averageZscoreOnsetIIE(delta,:),1,'omitnan');%bursts 26:34
%Theta:
averageOnsetT=mean(averageOnset(theta,:),1,'omitnan');
averageOnsetSeizureT=mean(averageOnsetSeizure(theta,:),1,'omitnan');%bursts 26:34
averageOnsetIIET=mean(averageOnsetIIE(theta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureT=mean(averageZscoreOnsetSeizure(theta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIET=mean(averageZscoreOnsetIIE(theta,:),1,'omitnan');%bursts 26:34
%
%Median over frequency bands:
%G1
averageOnsetG1=median(averageOnset(gamma1,:),1,'omitnan');
averageOnsetSeizureG1=median(averageOnsetSeizure(gamma1,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG1=median(averageOnsetIIE(gamma1,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG1=median(averageZscoreOnsetSeizure(gamma1,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG1=median(averageZscoreOnsetIIE(gamma1,:),1,'omitnan');%bursts 26:34
%G1_5
averageOnsetG1_5=median(averageOnset(gamma1_5,:),1,'omitnan');
averageOnsetSeizureG1_5=median(averageOnsetSeizure(gamma1_5,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG1_5=median(averageOnsetIIE(gamma1_5,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG1_5=median(averageZscoreOnsetSeizure(gamma1_5,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG1_5=median(averageZscoreOnsetIIE(gamma1_5,:),1,'omitnan');%bursts 26:34
%G2
averageOnsetG2=median(averageOnset(gamma2,:),1,'omitnan');
averageOnsetSeizureG2=median(averageOnsetSeizure(gamma2,:),1,'omitnan');%bursts 26:34
averageOnsetIIEG2=median(averageOnsetIIE(gamma2,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureG2=median(averageZscoreOnsetSeizure(gamma2,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEG2=median(averageZscoreOnsetIIE(gamma2,:),1,'omitnan');%bursts 26:34
%Alpha:
averageOnsetA=median(averageOnset(alpha,:),1,'omitnan');
averageOnsetSeizureA=median(averageOnsetSeizure(alpha,:),1,'omitnan');%bursts 26:34
averageOnsetIIEA=median(averageOnsetIIE(alpha,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureA=median(averageZscoreOnsetSeizure(alpha,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEA=median(averageZscoreOnsetIIE(alpha,:),1,'omitnan');%bursts 26:34
%Beta:
averageOnsetB=median(averageOnset(beta,:),1,'omitnan');
averageOnsetSeizureB=median(averageOnsetSeizure(beta,:),1,'omitnan');%bursts 26:34
averageOnsetIIEB=median(averageOnsetIIE(beta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureB=median(averageZscoreOnsetSeizure(beta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIEB=median(averageZscoreOnsetIIE(beta,:),1,'omitnan');%bursts 26:34
%Delta:
averageOnsetD=median(averageOnset(delta,:),1,'omitnan');
averageOnsetSeizureD=median(averageOnsetSeizure(delta,:),1,'omitnan');%bursts 26:34
averageOnsetIIED=median(averageOnsetIIE(delta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureD=median(averageZscoreOnsetSeizure(delta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIED=median(averageZscoreOnsetIIE(delta,:),1,'omitnan');%bursts 26:34
%Theta:
averageOnsetT=median(averageOnset(theta,:),1,'omitnan');
averageOnsetSeizureT=median(averageOnsetSeizure(theta,:),1,'omitnan');%bursts 26:34
averageOnsetIIET=median(averageOnsetIIE(theta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetSeizureT=median(averageZscoreOnsetSeizure(theta,:),1,'omitnan');%bursts 26:34
averageZscoreOnsetIIET=median(averageZscoreOnsetIIE(theta,:),1,'omitnan');%bursts 26:34
figure('DefaultAxesFontSize',14)
     hold on
     imagesc([-19:20],interp_frqs,averageZscoreOnsetIIE)%of onset/offset
     %imagesc([-19:20],interp_frqs,averageZscoreOnsetSeizure)%of onset/offset
     %caxis([-1.4385    1.0482])%seizures
     ylabel('Frequency (Hz)')
%      yticks([0 1 2])
%      yticklabels({'1','10','100'})
     %ylim([0 1.9])
     set(gca,'YDir','normal') 
     colorbar%'PSD (dB/Hz)'
     shading interp; % interpolate shading to make color transitions visually smooth
     title(sprintf('IIE median z score PSD onset'))
     title(sprintf('Seizures median z score PSD onset'))
     axis tight
xlabel('Time (sec)')
 xticks([-10 0 10])
     xticklabels({'-5','0','5'})
     caxis([-0.568 1.112])
      %xlim([-17 17])%-19.5 20.5
      saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Onsets/medianZscoreSeiz_Onset_burst%d.png',mouse_id,ip))
      %/bursts_PSDonsetIIE.png',ip));
      saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Onsets/medianZscoreIIE_Onset_burst%d.png',mouse_id,ip))
%
data1_mean = mean(averageZscoreOnsetIIE(delta,:),1,'omitnan');%67 freqs% Mean Across Columns
data2_mean = mean(averageZscoreOnsetSeizure(delta,:),1,'omitnan');%67 freqs% Mean Across Columns
data1_SEM = std(averageZscoreOnsetIIE(delta,:),1,'omitnan')/sqrt(length(delta));% SEM Across Columns
data2_SEM = std(averageZscoreOnsetSeizure(delta,:),1,'omitnan')/sqrt(length(delta));% S

data1_mean = mean(averageZscoreOnsetIIE(gamma1,:),1,'omitnan');%67 freqs% Mean Across Columns
data2_mean = mean(averageZscoreOnsetSeizure(gamma1,:),1,'omitnan');%67 freqs% Mean Across Columns
data1_SEM = std(averageZscoreOnsetIIE(gamma1,:),1,'omitnan')/sqrt(length(gamma1));% SEM Across Columns
data2_SEM = std(averageZscoreOnsetSeizure(gamma1,:),1,'omitnan')/sqrt(length(gamma1));% S

x1=data1_mean
x2=data2_mean
[h,p,ks2]=kstest2(x1,x2)
figure('DefaultAxesFontSize',14)
hold on
%   h1=subplot(2,1,1)% h1 = nexttile;
%   hold on
shadedErrorBar([-17:17], data1_mean, data1_SEM,'lineprops', '-b')
shadedErrorBar([-17:17], data2_mean, data2_SEM,'lineprops', '-r')
title('Delta (1-4 Hz) PSD onsets')
title('Gamma low (30-60 Hz) PSD onsets')
ylabel('PSD (dB/Hz)')
legend('IIEs','Seizures')
xlabel('Time (sec)')
axis tight
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Onsets/gamm1_Onsets.png',mouse_id))
saveas(gcf,sprintf('E:/figures/mouse%d/4AP/EEG/Buffers/Onsets/delta_Onsets.png',mouse_id))

%offsets figures loading PSD results:
for iseg=1:nseg 
    tstartEEG0=1+(iseg-1)*seg_frames;%EEG counter  time started from zero
    tendEEG0=iseg*seg_frames%
    tstartCa0=1+(iseg-1)*seg_framesCa;%EEG counter  time started from zero
    tendCa0=iseg*seg_framesCa%yinitial1=eegData.CleanEEG_cV(tstart0:tend0, 1);%whole segment   
    %lips=find(xxx(:,1)==iseg)%1:6/7,8 start at first segment/2 at second
    lips=find(xxx(:,2)==iseg)% %around onset and around offset:
   % lips=lips2;
    for ipp=1:length(lips)%1:length(plateau_start)%segments 8,10 without bursts
        ip=lips(ipp);%ip2=lips2(ipp);
        tstart=round(((iseg-1)*(-CaSeg+EEGseg)+plateau_end(ip))*337);%=round(((iseg-1)*(-CaSeg+EEGseg)+plateau_end(ip))*337);    
        tstartEEG=tstart-20*bin+1;
        tstartCa=plateau_end(ip)-round(20*binCa+1);
        if tstartCa<0 %then also tstartEEG<0
            tstartCa=tstartCa0;%segment start
            tstartEEG=tstartEEG0;
        end    
        t_endEEG=tstart+20*bin;
        t_endCa=plateau_end(ip)+round(20*binCa);
        %symmetric to onset:
        if t_endCa > tendCa0
            t_endCa=tendCa0;%segment end
            t_endEEG=tendEEG0;
        end
         load(sprintf('../../Results/Buffers/OnsetsOffsets/Null_4AP/meanRandom_%dSeg_matrixOffset_burst%d.mat',iseg,ip))
        %;,'meanRandVmatrixOnset');
        load(sprintf('../../Results/Buffers/OnsetsOffsets/Null_4AP/stdRandom_%dSeg_matrixOffset_burst%d.mat',iseg,ip));
        %,'stdRandVmatrixOnset');
        load(sprintf('../../Results/Buffers/OnsetsOffsets/Null_4AP/semRandom_%dSeg_matrixOffset_burst%d.mat',iseg,ip));
        %,'semRandVmatrixOnset');
        load(sprintf('../../Results/Buffers/OnsetsOffsets/Null_4AP/zscore_%dSeg_matrixOffset_burst%d.mat',iseg,ip));
        %,%z_scoreV=(matrixOnset-meanRandVmatrixOnset)./stdRandVmatrixOnset;
        figure('DefaultAxesFontSize',14)
        hold on
        %find the exact offset location and plot it in bold:
        arx=floor(2*(plateau_end(ip)-tstartCa)*0.0337 -2.5);%10
        tel=floor(2*(t_endCa-plateau_end(ip))*0.0337-2.5);%17
        imagesc([-arx:tel].*0.5,interp_frqs,meanRandmatrixOffset);%
        ylabel('Frequency (Hz)')
     yticks([0 1 2])
     yticklabels({'1','10','100'})
     caxis% 0.1998   117.8
     shading interp; % interpolate shading to make color transitions visually smooth
     title(sprintf('Burst %d offset random mean (4AP null)',ip))
      xlabel('Time (sec)')
      axis tight
      xlim([-10 10])
      colorbar
      saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/meanRandom_%dSeg_matrixOffset_burst%d.png',iseg,ip));
        saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/meanRandom_%dSeg_matrixOffset_burst%d.fig',iseg,ip));
        close
        figure('DefaultAxesFontSize',14)
        hold on
        imagesc([-arx:tel].*0.5,interp_frqs,semRandmatrixOffset);%
        ylabel('Frequency (Hz)')
     yticks([0 1 2])
     yticklabels({'1','10','100'})
     caxis% 0.1998   117.8
     shading interp; % interpolate shading to make color transitions visually smooth
     title(sprintf('Burst %d offset SEM (4AP null)',ip))
      xlabel('Time (sec)')
      axis tight
      xlim([-10 10])
      colorbar
       saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/semRandom_%dSeg_matrixOffset_burst%d.png',iseg,ip));
        saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/semRandom_%dSeg_matrixOffset_burst%d.fig',iseg,ip));
        close
        figure('DefaultAxesFontSize',14)
        hold on
        imagesc([-arx:tel].*0.5,interp_frqs,stdRandmatrixOffset);%
        ylabel('Frequency (Hz)')
     yticks([0 1 2])
     yticklabels({'1','10','100'})
     caxis% 0.1998   117.8
     shading interp; % interpolate shading to make color transitions visually smooth
     title(sprintf('Burst %d offset STD (4AP null)',ip))
      xlabel('Time (sec)')
      axis tight
      xlim([-10 10])
      colorbar
       saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/stdRandom_%dSeg_matrixOffset_burst%d.png',iseg,ip));
        saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/stdRandom_%dSeg_matrixOffset_burst%d.fig',iseg,ip));
        close
        figure('DefaultAxesFontSize',14)
        hold on
        imagesc([-arx:tel].*0.5,interp_frqs,z_score);%
        ylabel('Frequency (Hz)')
      yticks([0 1 2])
     yticklabels({'1','10','100'})
     caxis% 0.1998   117.8
     shading interp; % interpolate shading to make color transitions visually smooth
     title(sprintf('Burst %d offset z score (4AP null)',ip))
      xlabel('Time (sec)')
      axis tight
      xlim([-10 10])
      colorbar
       saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/zscore_%dSeg_matrixOffset_burst%d.png',iseg,ip));
        saveas(gcf,sprintf('../../SingleImages/Buffers/OnsetsOffsets/Null_4AP/zscore_%dSeg_matrixOffset_burst%d.fig',iseg,ip));
        close
    end
end


avrgOnset=nan(nfreqs,35,100);
avrgMean=nan(nfreqs,35,100);
avrgSTD=nan(nfreqs,35,100);
avrgSEM=nan(nfreqs,35,100);
avrgZscore=nan(nfreqs,35,100);
