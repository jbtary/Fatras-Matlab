%% Real data given by JB.

clear
addpath(genpath('crewes'));
load ('rawdata.mat');
%% Plot data.

SegyData=SegyData(:,78:-1:1);

%%
plotimage(SegyData',time,depth);hold on;xlabel('Time (ms)'); ylabel('Depth (m)');title('Total Downgoing Field');

%% First arrival times.
for i=1:length(SegyData(1,:))
    
    index = find(SegyData(:,i)~=0, 1, 'first');
    tpp=time(index);
    ind(i,1)=index;
    tp(i,1)=tpp;
end
%%
%% Make dt array between receivers. first.time.pick(i+1)-firts.time.pickt(i).
% ddt is the travel time difference for the firts time picks between receivers.
for i=1:length(tp)-1
d=tp(i+1)-tp(i);
ddt(:,i)=d;end
%% Plot time.

%depth=depth';
depth=depth(:,78:-1:1);

%%
figure;
plot(depth,tp,'o');hold on; plot(depth,tp_bpicks,'o');xlabel('Time (ms)');ylabel('Depth (m)'); title('Arrival time for each receiver'); legend('Tary et al. (2017)','Estimated time')
%% Separate.
[vspdown,vspup,vspupflat,tpup,rdown]=VSP_separation(SegyData,time,depth,tp_bpicks,6,6);
%[vspdown,vspup,vspupflat,tpup,rdown]=VSP_separation(vspdown,time,depth,tp_bpicks,6,6);
%%
plotimage(vspdown');hold on;xlabel('Time (ms)'); ylabel('Depth (m)');title('Total Primary Downgoing Field');
%plotimage(vspup');xlabel('Time (ms)'); ylabel('Receiver #');title('Upgoing Field');
%% Plot Original waves.
figure; hold on;
for i=1:length(SegyData(1,:))
plot(time,SegyData(:,i)); end
%%
figure;
plot(time,SegyData(:,8));hold on; plot(time,vspdown(:,8),'red');
%% Plot downgoinf field.
figure; hold on;
for i=1:5:length(vspdown(1,:))
plot(time,vspdown(:,i)); 
end
%%
figure;
plot(time,vspdown(:,50),'black');hold on;
plot(time,SegyData(:,50));hold on;
vspdown(ind(50)+100:end,50)=0; 
plot(time,vspdown(:,50));
%% Plot frecuency.
nfft=1024;
x=(1:nfft/2)*(fs/nfft);
figure; hold on;xlabel('Frequency (Hz)');ylabel('FFT'); title('Frequency espectro for each receiver')
for i=1:length(vspdown(1,:))
     magfreq=fft(vspdown(1:end,i),nfft);
     magfreq=abs(magfreq(1:end/2));
     plot(x,magfreq);
end
%% 
figure;
    magfreq=fft(vspdown(ind(8):ind(8)+56,8),nfft);
     magfreq=abs(magfreq(1:end/2));
     plot(x,magfreq);
%%

for i=1:78
    vspdown(ind(i)+100:end,i)=0;    
end
%% Filter Data
fNy=fs/2;
[b,a]=butter(6,[5 45]/fNy,'bandpass'); % Butterworth filter of order 6
for i=1:length(vspdown(1,:))    
filterdata(:,i)= filtfilt(b,a,vspdown(:,i));
end


%% Plot filteres freq.

figure; hold on;
for i=1:length(vspdown(1,:))
     magfreq=fft(filterdata(ind(i):end,i),nfft);
     magfreq=abs(magfreq(1:end/2));
     plot(x,magfreq);
end

%%

[QSR,QSRtrace,yy,x2,freqq]=QSpectralRadios(vspdown,fs,tp_bpicks,15,50,25,26,depth,0,ind)
%%

[QSRD,QSRtraceD,yyD,x2D,freqqD]=DCIQSRSpectralRadios(vspdown,fs,tp_bpicks-tp_bpicks(1),15,50,25,26,depth,5,0,ind)
%%
figure;
plot(1./QSR, depth (1:length(QSR)));flipy; hold on;
plot(1./QSRD, depth (1:length(QSRD)));
%%
figure;
plot(QSR); hold on; plot(QSRD);
%%
% MMFS method.(Knowing a real Q value of index 4 depth. (Absolute attenuation).
% and assuming an hypothetical Q value of index 4 depth. (Relative attenuations).
%tteo=tteo(1:132);
[QMMFS,QMMFS1,step4]=MMFS(yyD,ddt(2:76),x2D,4,60,100,depth);
%V=errperf(Q(2:length(QMMFS)),QMMFS(2:length(QMMFS)),'rmspe')  % Root Mean Squared Error
%% Find scaterring attenuation
[QMMFS2, Intri]=MMFS2(yyD,ddt(2:76),tp-(tp(1)),x2D,depth,freqqD,7);
figure;
plot((1./QMMFS2(1:(length(QMMFS2)))),depth(1:(length(QMMFS2))));flipy; xlabel('1/Q'); ylabel('depth (m)'); title('Absolute QMMFS estimation');

%%
[f,t,Intrinsecf,Qestt,A,QMMFSI,TOTAL]=MMFS3(yyD,ddt(2:76),tp-(tp(1)),x2D,depth,freqqD,7)%,dp1,dp2,dp3,dp4)

 %%
 figure; hold on;flipy;xlabel('1/Qs'); ylabel('Depth'); title('Qs values vs depth for each frequency')
 for i=1:36
     plot(1./QMMFSI(i,3:end),depth(3:length(QMMFSI)));
     
 end
