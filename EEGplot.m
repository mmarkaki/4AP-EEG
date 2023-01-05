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

avrgOnset=nan(nfreqs,35,100);
avrgMean=nan(nfreqs,35,100);
avrgSTD=nan(nfreqs,35,100);
avrgSEM=nan(nfreqs,35,100);
avrgZscore=nan(nfreqs,35,100);
