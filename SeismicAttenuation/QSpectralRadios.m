function[QSR,QSRtrace,yy,x2,freqq]=QSRSpectralRadios(data,fs,tp,mf,hf,trace1,trace2,z,pp,ind,Q)

nfft=1024;    
x=fs*(1:(nfft/2))/nfft;
x2=x(1,mf:hf);
    
for i=1:(length(data(1,:))-2)
    
    hhf=hf;
    %if (i<15)
     %   hhf=200;
    %end
    magfreQSR=0;
    magfreQSR2=0;
    
    magfreQSR=fft(data(ind(i):end,i),nfft);
    magfreQSR=2*abs(magfreQSR(1:end/2));
    freq=log(magfreQSR(mf:hhf,1)./2);
    
    magfreQSR2=fft(data(ind(i+1):end,i+1),nfft);
    magfreQSR2=2*abs(magfreQSR2(1:end/2));
    
    y=log(magfreQSR2./magfreQSR);
    y2=y(mf:hhf,1);
    
    p = polyfit(x2,y2',1);
    f = polyval(p,x2);
    
    %Find the final atenuation.
    QSRest=(-pi*(tp(i+1)-tp(i)))/p(1);
    
    
    listmagfrec(:,i)=magfreQSR;
    listmagfrec2(:,i)=magfreQSR2;
    yy(:,i)=y2;
    freqq(:,i)=freq;
    QSR(:,i)=abs(QSRest);
    
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
subplot(3,1,1); plot(x,listmagfrec(:,trace1));ylabel('FFT 1');title('Slope found in the SR method');
subplot(3,1,2); plot(x,listmagfrec2(:,trace2));ylabel('FFT 2');
subplot(3,1,3); plot(x,listy,listx2,listf,'-');ylabel('log(FFT2/FFT1)');xlabel('Frecuency (Hz)');
%QSR=QSR(:,length(QSR):-1:1);

%V=errperf(Q(3:length(QSR)),QSR(3:length(QSR))','rmspe')  % Root Mean Squared Error
figure;
plot((1./QSR(1:(length(QSR)-pp))),z(1:(length(QSR)-pp)));flipy; xlabel('1/Q'); ylabel('depth (m)'); title('QSR estimation');
hold on;plot(1./Q(1:length(Q)),z(1:length(Q)));legend('Estimathed Q by SR.','Real Q.')


