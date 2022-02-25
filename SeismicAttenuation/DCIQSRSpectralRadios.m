function[QSR,QSRtrace,yy,x2,freqq]=DCIQSRSpectralRadios(data,fs,tp,mf,hf,trace1,trace2,z,Q,pp,ind)

nfft=1024;    
x=fs*(1:(nfft/2))/nfft;
x2=x(1,mf:hf);
    magfreQSR0=fft(data(ind(1):end,1),nfft);
    magfreQSR0=abs(magfreQSR0(1:end/2));
eps=0.1;
    
for i=2:(length(data(1,:))-2)
    
    hhf=hf;
    %if (i<15)
     %   hhf=200;
    %end
    magfreQSR=0;
    magfreQSR2=0;
    
    magfreQSR=fft(data(ind(i):end,i),nfft);
    magfreQSR=2*abs(magfreQSR(1:end/2));
    %magfreQSR=(magfreQSR)./magfreQSR0;
    %magfreQSR=(magfreQSR.*conj(magfreQSR0))./(magfreQSR0.*conj(magfreQSR0)+eps*mean(magfreQSR0.*conj(magfreQSR0)));
    magfreQSR=(magfreQSR./magfreQSR0);%FOR DECONVOLUTION
    % magfreQSR=(magfreQSR.*conj(magfreQSR0)); % FOR CROSSCORRELATION
    
    freq=log(magfreQSR(mf:hhf,1)./2);
    
    magfreQSR2=fft(data(ind(i+1):end,i+1),nfft);
    magfreQSR2=2*abs(magfreQSR2(1:end/2));
    %magfreQSR2=magfreQSR2./magfreQSR0;
    %magfreQSR2=(magfreQSR2.*conj(magfreQSR0))./(magfreQSR0.*conj(magfreQSR0)+eps*mean(magfreQSR0.*conj(magfreQSR0)));
    magfreQSR2=(magfreQSR2./magfreQSR0); %FOR DECONVOLUTION
    %magfreQSR2=(magfreQSR2.*conj(magfreQSR0)); %FOR CROSSCORRELATION
    
    y=log(magfreQSR2./magfreQSR);
    y2=y(mf:hhf,1);
    
    p = polyfit(x2,y2',1);
    f = polyval(p,x2);
    
    QSRest=(-pi*(tp(i+1)-tp(i)))/p(1);
    
    
    listmagfrec(:,i)=magfreQSR;
    listmagfrec2(:,i)=magfreQSR2;
    yy(:,i-1)=y2;
    freqq(:,i-1)=freq;
    QSR(:,i-1)=abs(QSRest);
    
end 
%if (trace2 < 15)
%    hf=200;
%end
listy=log(listmagfrec(:,trace2)./listmagfrec(:,trace1));
listy2=listy(mf:hf,1);
listx2=x(1,mf:hf);
listp = polyfit(listx2,listy2',1);
listf = polyval(listp,listx2);
QSRtrace=(-pi*(tp(trace2)-tp(trace1)))/listp(1);
disp(QSRtrace);
disp(listp(1));
figure;  
subplot(3,1,1); plot(x,listmagfrec(:,trace1));ylabel('FFT 1');;title('Slope found in the SR method');
subplot(3,1,2); plot(x,listmagfrec2(:,trace2));ylabel('FFT 2');
subplot(3,1,3); plot(x,listy,listx2,listf,'-');ylabel('log(FFT2/FFT1)');xlabel('Frecuency (Hz)');
%QSR=QSR(:,length(QSR):-1:1);

%V=errperf(Q(3:length(QSR)),QSR(3:length(QSR))','rmspe')  % Root Mean Squared Error
figure;
plot((1./QSR(1:(length(QSR)-pp))),z(1:(length(QSR)-pp)));flipy; xlabel('1/Q'); ylabel('depth (m)'); title('DCI+QSR estimation');
%hold on;plot(1./Q(1:length(Q)),z(1:length(Q)));legend('Estimathed Q by SR.','Real Q.')