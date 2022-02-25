%%PlOT data

plotimage(downq',tq,zr);hold on; 
xlabel('Time (ms)'); ylabel('Receiver Depth (m)');title('Donwoing Field');

%% Theoretical time arrival.

teo=0;
   for i=1:length(zr)-1
       
     m=1500;
   if i> 19;
       m=1000;
   end
   
   if i> 74;
      m=2000;
   end
   
   if i> 133;
       m=1000;
   end
   
    teo=teo+((zr(i+1)-zr(i))/m);
   
   tteo(:,i)=teo;
   end

%% Make dt array between receivers. first.time.pick(i+1)-firts.time.pickt(i).
% ddt is the travel time difference for the firts time picks between receivers.
for i=1:length(tteo)-1
d=tteo(i+1)-tteo(i);
ddt(:,i)=d;
end
%% Find times.(Finite difference method). ind is the index of the first pick. 
%Note that the index can not be calculated theoricaly because the tteo values does not match with any specific value of the time vector.
% This index (ind) is very importat because the FFT is going to be taken from this to the end of the array.

[tp,ind]=PickTime(downq,tq,dt);
%% Plot time.

figure;
plot(z(1:length(tteo)-1),tteo(1:length(tteo)-1),'o');hold on;
plot(z(1:length(tp)-1),tp(1:length(tp)-1),'o');legend('Theoric time','Estimated time');
xlabel('depth (m)');ylabel('First arrival time (s)');title('First arrival time for each receiver at depth z');
%% Create theorical Q list.
Q=zeros(1,length(zr));
Q(1:19)=30;
Q(19:74)=15;
Q(74:133)=30;
Q(133:length(zr))=40;
%% PLot waves.

figure; hold on; xlabel('Time (s)'); title('Waves');
for i=40:2:50
plot(tq,vspq(:,i));
end
%%
figure;
plot(tq,vspq(:,60));hold on;
plot(tq,downq(:,60));hold on;
downq(ind(60)+84:end,60)=0;
plot(tq,downq(:,60));
%%

for i=1:149
    downq(ind(i)+84:end,i)=0;    
end
%%
figure;
plot(tq,downq(:,130));
%%
figure; hold on;
for i=1:149
    plot(downq(:,i));    
end
    %%
h=10;
figure;
plot(downq(:,h));
%% Check frecuency spectrum.
% From the arrival time index to the end of the array. 
fs=1/dt;
nfft=1024;
x=(1:nfft/2)*(fs/nfft);
figure; hold on; xlabel('Frecuency (Hz)');ylabel('FFT'); title('Frecuency espectrum for each receiver at depth');
for i=19:5:(length(downq(1,:)))
    
    magfreq=abs(fft(downq(ind(i):end,i),nfft));
    magfreq= magfreq(1:end/2);
    
    plot(x(1:180),magfreq(1:180));
end 
%% Separate.
[vspdown,vspup,vspupflat,tpup,rdown]=VSP_separation(vspq(:,1:149),tq,z,tp,6,6);
%[vspdown2,vspup2,vspupflat,tpup,rdown]=VSP_separation(vspdown,tq,z,tteo,6,6);

 plotimage(vspup2');

 plotimage(vspdown2');
%% Filter the data.

fNy=fs/2;
[b,a]=butter(6,[1 100]/fNy,'bandpass'); % Butterworth filter of order 6
for i=1:length(downq(1,:))    
filterdata(:,i)= filtfilt(b,a,downq(:,i));
end

%% Q by SR method.
l1=downq(:,19:end);
l2=tteo(19:end);
l3=zr(19:end);
l4=Q(19:end);
l5=ind(19:end);
[QSR,QSRtrace,yy,x2,freq]=QSpectralRadios(downq,1/dt,tteo,15,60,135,140,zr,Q,0,ind,36);
V=errperf(Q(20:length(QSR)),QSR(20:length(QSR)),'rmspe')  % Root Mean Squared Error
%%
mf=15;
hf=60;
nfft=1024;    
x=fs*(1:(nfft/2))/nfft;
x2=x(1,mf:hf);
    
for i=1:(length(l1(1,:))-2)
    
    hhf=hf;
    %if (i<15)
     %   hhf=200;
    %end
    magfreQSR=0;
    magfreQSR2=0;
    
    magfreQSR=fft(l1(l5(i):end,i),nfft);
    magfreQSR=2*abs(magfreQSR(1:end/2));
    freq=log(magfreQSR(mf:hhf,1)./2);
    
    magfreQSR2=fft(l1(l5(i+1):end,i+1),nfft);
    magfreQSR2=2*abs(magfreQSR2(1:end/2));
    
    y=log(magfreQSR2./magfreQSR);
    y2=y(mf:hhf,1);
    
    p = polyfit(x2,y2',1);
    f = polyval(p,x2);
    
    QSRest=(-pi*(l2(i+1)-l2(i)))/p(1);
    
    
    listmagfrec(:,i)=magfreQSR;
    listmagfrec2(:,i)=magfreQSR2;
    yy(:,i)=y2;
    freqq(:,i)=freq;
    QSRR(:,i)=abs(QSRest);
    
end
%%
figure;
plot(QSR);
%%
figure;
plot((1./QSRR(1:(length(QSRR)))),l3(1:(length(QSRR))),'o');flipy; xlabel('1/Q'); ylabel('depth (m)'); title('QSR estimation');
hold on;plot(1./l4(1:length(QSRR)),l3(1:length(QSRR)));legend('Estimathed Q by SR.','Real Q.')
%% MMFS method.

[QMMFS,QMMFS1,step4,A]=MMFS(yy,ddt,x2,4,30,5,150,zr,Q);
V=errperf(Q(2:length(QMMFS)),QMMFS(2:length(QMMFS)),'rmspe')  % Root Mean Squared Error

%% MMFS method to find absolute attenuations without asuming any hypothetical Q value.

[QMMFS2, Intri,A2]=MMFS2(yy,ddt,tteo,x2,zr,freq,Q);
P0=1./QMMFS2;
P2=1./Q;
figure;
plot((1./QMMFS2(1:(length(QMMFS2)))),zr(1:(length(QMMFS2))),'o');flipy; xlabel('1/Q'); ylabel('depth (m)'); title('Absolute QMMFS estimation');
hold on;plot(1./Q(1:length(QMMFS2)),zr(1:length(QMMFS2)));
;legend('Absulute Q by MMFS.','Real Q.','Relative Q by MMFS.')
V=errperf(Q(2:length(QMMFS2)),QMMFS2(2:length(QMMFS2)),'rmspe')  % Root Mean Squared Error