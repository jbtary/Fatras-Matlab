%% Syntetick data given by JB. (RUN IT IN ORDER).

clear
addpath(genpath('crewes'));
load ('synth_VSP2.mat')

%% Plot data.
datac=datac';
plotimage(datac,zarr);xlabel('Time (ms)'); ylabel('Depth (m)'), title('Synthetick field');

%% plot waves.

figure; hold on; xlabel('Time (s)'); title('Waves');
for i=1:4:length(zarr)
plot(times(1:1000),datac(1:1000,i));end
%%
figure; plot(times(1:1000),datac(1:1000,40));
%% Theorical firts pick times.

teo=0;
   for i=1:length(zarr)-1
       
     m=1500;
   if i> 19
       m=1000;
   end
   
   if i> 74
      m=2000;
   end
   
   if i> 119
       m=1000;
   end
   
    teo=teo+((zarr(i+1)-zarr(i))/m);
   
   tteo(:,i)=teo;end
%% Make dt array between receivers. first.time.pick(i+1)-firts.time.pickt(i).
% ddt is the travel time difference for the firts time picks between receivers.
for i=1:length(tteo)-1
d=tteo(i+1)-tteo(i);
ddt(:,i)=d;end

%% Find times.(Finite difference method). ind is the index of the first pick. 
%Note that the index can not be calculated theoricaly because the tteo values does not match with any specific value of the time vector.
% This index (ind) is very importat because the FFT is going to be taken from this to the end of the array.

[tp,ind]=PickTime(datac,times,dt);
%% Theoretical index of first pick. THIS WONT WORK.
for i=1:length(tteo(1,:))
    
    index = find(times(1,:)==tteo(i));
    indt(i,1)=index;
end

%% Plot time.
% Note that the times are very similar, so the pick time method is well done.
figure;
plot(zarr(1:length(zarr)-1),tteo(1:length(zarr)-1),'o');hold on;
plot(zarr(1:length(zarr)-1),tp(1:length(zarr)-1),'o');legend('Theoric time','Estimated time');
xlabel('depth (m)');ylabel('First arrival time (s)');title('First arrival time for each receiver at depth z');
%% Create theorical Q list.
Q=zeros(1,length(zarr));
Q(1:19)=60;
Q(19:74)=30;
Q(74:119)=60;
Q(119:length(zarr))=80;
%% Check frecuency spectrum.
% From the arrival time index to the end of the array. 
nfft=1024;
x=(1:nfft/2)*(fs/nfft);
figure; hold on; xlabel('Frecuency (Hz)');ylabel('FFT'); title('Frecuency espectrum for each receiver at depth');
for i=1:(length(datac(1,:)))
    
    magfreq=abs(fft(datac(ind(i):end,i),nfft));
    magfreq= magfreq(1:end/2);
    
    plot(x(1:180),magfreq(1:180));
end 

%% Q by SR method.
[QSR,QSRtrace,yy,x2,freq]=QSpectralRadios(datac,1000,tteo,20,80,1,2,zarr,0,ind,Q);
V=errperf(Q(2:length(QSR)),QSR(2:length(QSR)),'rmspe')  % Root Mean Squared Error

%% Check constant (B) with the Q values given by the SR method.
% With the graph is very clear that in a given medium with constant Q, the
% factor B depend very vey litle on depth ( almost an straight line for a
% constant Q), but when there is a change in the medium with a different Q 
%the position of the line changes completely. This mean that the B factor
%depends weakly on depth and strongly on the frecuency content (Q changes).
for i=1:length(ddt)
c=exp(mean(yy(:,i)+x2'.*pi*ddt(i)./QSR(i)));
cd=exp(mean(yy(:,i)+x2'.*pi*ddt(i)./Q(i)));
cc(:,i)=c;
ccd(:,i)=cd;
end
figure; 
plot(cc(2:end),zarr(2:length(cc)));flipy; xlabel('Factor (B)');ylabel('depth (m)'); title('mean Factor (B) vs Depth');
hold on;plot(ccd(2:end),zarr(2:length(ccd))); legend('With estimated QSR','With real Q');
%% Plot factor B over frequencies.
figure; hold on;xlabel('Frequency (Hz)'); ylabel('A factor'); title('A factor vs Frequency for each receiver at depth');
for i=1:131
c=exp((yy(:,i)+x2'.*pi*ddt(i)./QSR(i)));
cd=exp((yy(:,i)+x2'.*pi*ddt(i)./Q(i)));
ccf(:,i)=c;
ccdf(:,i)=cd;
plot(x2,ccdf(:,i+1));end
%% ploting Factor B vs depth not for the mean factor B but for all the frequencies.
figure; hold on;xlabel('Depth (m)');ylabel('B factor');title('B factor vs depth for each frequency');
for i=1:51
plot(zarr(2:132),ccdf(i,2:end));end

%% Add noise and SR method.
% Without filtering the data we can note that the SR method is very
% sensitive to noise. 
[noisydata]=addnoisy(datac,25,1);

figure; 

plot(times,noisydata(:,40));hold on;
plot(times,datac(:,40))


[QSR,QSRtrace,yy,x2,freqq]=QSpectralRadios(noisydata,1000,tteo,10,60,50,53,zarr,0,ind,Q);
V=errperf(Q(2:length(QSR)),QSR(2:length(QSR)),'rmspe')  % Root Mean Squared Error
%% SR method plus DC interferometry.

[fdata]=InterferometryDC(datac,1000,15,100,ind);
figure;
plot(fdata(:,1));
%%
tteo2=tteo-tteo(1);
%%
for i=1:length(tteo2)-1
d=tteo2(i+1)-tteo2(i);
ddt2(:,i)=d;end
%%
[tp2,ind2]=PickTime(fdata,times,dt);
figure;
plot(zarr(1:length(zarr)-1),tp2(1:length(zarr)-1));
%%
figure; hold on
for i=1:length(zarr)-1
plot(fdata(:,i)); end

%%
[QSRS,QSRtraceS,yyS,x2S,freqS]=DCIQSRSpectralRadios(datac,1000,tteo,10,60,30,35,zarr,Q,0,ind);
V=errperf(Q(5:length(QSRS)),QSRS(5:length(QSRS)),'rmspe')  % Root Mean Squared Error

%% SR method plus CC interferometry.

[xcdfinal,lagfinal,indb]=InterferometryCC(datac,1000,15,80,50);
[QSR,QSRtrace,RMSE]=QSpectralRadios(xcdfinal,1000,tteo,10,60,50,53,zarr,Q,50);
%% MMFS method.(Knowing a real Q value of index 4 depth. (Absolute attenuation).
% and assuming an hypothetical Q value of index 4 depth. (Relative attenuations).
tteo=tteo(1:132);
[QMMFS,QMMFS1,step4,step3,bet,step2,step1]=MMFS(yy,ddt,x2,4,60,30,zarr,Q);
V=errperf(Q(2:length(QMMFS)),QMMFS(2:length(QMMFS)),'rmspe')  % Root Mean Squared Error

%% Plot factor A resulting with the MMFS method using the real Q.
% Note that this is less than the B given by the SR method. The dependency on frecuency is much less.
step11=mean(bet);
cc2=step4+(ddt./Q(1:length(ddt)));
cc3=step4+(ddt./QMMFS(1:length(ddt)));

%%
figure; 
plot(cc2(2:end),zarr(2:length(cc3)));flipy; xlabel('Factor (A)');ylabel('depth (m)'); title('Factor (A) vs Depth');
title('Factor A given by MMFS method'); 
%% ploting Factor A vs depth the for the total phi portion get with step 3 (without making the last median frequency step)
figure; hold on;ylabel('Depth (m)');xlabel('A factor');title('A factor vs depth for each frequency');
ccd=bet+(ddt./Q(1:length(ddt)));
%ccdf=(bet-step2)+(ddt./Q(1:length(ddt)));

for i=1:61
plot(ccd(i,2:end),zarr(2:132));flipy;
end
%% Find scaterring attenuation 
Scateoric=zeros(1,132);


bet=yy./(pi.*x2');

step1=mean(bet);

step2=(median((bet-step1)'))';

step3=bet-step2;

step4=median(step3);


sca=step3-step4;

figure; hold on;xlabel('1/Qs'); ylabel('Depth (m)');title('Scaterring attenuation Qs value vs depth')
QMMFSS=-ddt./(sca-ccd(:,4));
for i=1:61
    
plot(1./QMMFSS(i,2:end),zarr(2:132));flipy;
    
end

hold on; plot(Scateoric,zarr(1:132),'black');legend('Estimated Qs','Theoric Qs')
%%
V=errperf(Scateoric,QMMFS(1,:),'rmspe')  % Root Mean Squared Error
%%
a=[1e-30 1e-30 1e-30 1e-30]
b=[0 2 0 3]
V=errperf(a,b,'rmspe')  % Root Mean Squared Error
%%
%% MMFS method to find absolute attenuations without asuming any hypothetical Q value.

[QMMFS2, Intri]=MMFS2(yy,ddt(1:132),tteo,x2,zarr,freqq,Q);
figure;
plot((1./QMMFS2(1:(length(QMMFS2)))),zarr(1:(length(QMMFS2))),'o');flipy; xlabel('1/Q'); ylabel('depth (m)'); title('Absolute QMMFS estimation');
hold on;plot(1./Q(1:length(QMMFS2)),zarr(1:length(QMMFS2)));
;legend('Absulute Q by MMFS.','Real Q.','Relative Q by MMFS.')
V=errperf(Q(2:length(QMMFS2)),QMMFS2(2:length(QMMFS2)),'rmspe')  % Root Mean Squared Error
%%
figure;
plotimage(Intri);