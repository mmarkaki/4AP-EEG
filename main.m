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
