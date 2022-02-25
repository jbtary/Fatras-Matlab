function[fdata]=InterferometryDC(data,fs,fm,fh,ind)
% eps: Parametro de regularizacion propuesto.
eps=0.1;
nfft=1024;  
n=length(data(:,1));
fNy = fs/2;
for i=1:(length(data(1,:))-1)
    
W = 1;

magfreq= (fft(W.*data(ind(i):end,1),n));
%magfreq=2*abs(magfreq(1:end/2));
magfreq2=(fft(W.*data(ind(i+1):end,i+1),n));
%magfreq2=2*abs(magfreq2(1:end/2));


sfinal=ifft(magfreq2./magfreq);
%sfinal=real(ifft((magfreq2.*conj(magfreq))./((magfreq.^2)+eps)));

fdata(:,i)=sfinal;

end
