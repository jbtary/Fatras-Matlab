function[xcdfinal,lagfinal,indb]=InterferometryCC(data,fs,fcl,fch,n)

dt = 1/fs; fNy = fs/2;

norm1b = 1;

[b,a]=butter(6,[fcl fch]/fNy,'bandpass'); % Butterworth filter of order 6

funcle = 0.4; % Length in seconds of NCC
segle = funcle/dt;% length of data to NCC

noverlap = round(0.5*segle); % 50% overlap used

indb = 1:(segle-noverlap):length(data(:,1))-(segle-1);
inde = segle:(segle-noverlap):length(data(:,1));
kk=1;
mm=1;
for i=1:(length(data(1,:))-1)
    
    %data(:,i)= filtfilt(b,a,data(:,i));
    %data(:,i+1)=filtfilt(b,a,data(:,i+1));
    
    if norm1b == 1
            %         data1 = sign(data1);
            %         data2 = sign(data2);
            % Alternative: normalize by envelope (L. Stehly)
            envd1 = hilbert(data(:,1)); envd1 = abs(envd1);
            envd2 = hilbert(data(:,i+1)); envd2 = abs(envd2);
            data(:,1) = data(:,1)./envd1;
            data(:,i+1) = data(:,i+1)./envd2;
            clear envd1 envd2
    end
    ll=1;
    jj=1;
    for nn = 1:length(indb)    
    dat1=data(:,1);
    dat2=data(:,i+1);
    [xctmp,lag] = xcorr(dat2(indb(nn):inde(nn),1),dat1(indb(nn):inde(nn),1),'coeff');
 
    xctmp = xctmp/max(abs(xctmp));
    
    xcdtmp(:,ll) = xctmp;
            ll = ll+1;
            
    lag2(:,jj) = lag;
            jj = jj+1;
    end
    
xcd = mean(xctmp,2,'omitnan');
lag3=mean(lag2,2,'omitnan');

xcdfinal(:,kk)=xcd;   
xcdfinal(:,kk) = xcdfinal(:,kk)/max(abs(xcdfinal(:,kk)));
kk=kk+1;

lagfinal(:,mm)=lag3;   
mm=mm+1;
end 

figure;
 
ax = axes('Position',[0.1 0.7 0.8 0.2]);
plot(lagfinal(:,n)*(1/fs),xcdfinal(:,n),'k'); ylabel('Avg. coeff.');xlabel('Lag time (s)')
%plot(time,xcdfinal(:,75),'k'); ylabel('Avg. coeff.');xlabel('Lag time (s)')
set(gca,'XTickLabel',[])
linkaxes(ax,'x');
