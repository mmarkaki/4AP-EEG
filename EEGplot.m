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
Fc = 60;%to be filtered out from EEG
Wc = Fc/(Fs/2);
BW = Wc/10;%35;
[b,a] = iirnotch(Wc(1),BW(1));%the filter characteristics
eegTimeSegment = eeg4AP.Time_s(tstart0:tend0);
yinitial=eeg4AP.CleanEEG_cV(tstart0:tend0, 1);%eg, we initially consider a whole segment
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
        %fb = cwtfilterbank('SignalLength',length(y1),'SamplingFrequency',Fs,'FrequencyLimits',[minfreq maxfreq], 'Wavelet',"amor");
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
