function [QMMFS2, Intri,A]=MMFS2(yy,ddt,tteo,x,z,freq,Q)%,dp1,dp2,dp3,dp4)
%between receivers.
bet=yy./(pi.*x');

step1=mean(bet);

step2=(median((bet-step1)'))';

step3=bet-step2;

step4=median(step3);

%from source to recerivers. To calculate test Q.
betf=freq./(pi.*x');

step1f=mean(betf);

step2f=(median((betf-step1f)'))';

step3f=betf-step2f;

step4f=median(step3f);

%A=step4(dep)+(ddt(dep)./realQ);

for i=1:8
    f(:,i)=mean(step4f(3:i+3));
    t(:,i)=mean(tteo(3:i+3));
end 
p=polyfit(t,f,1);
Qest=-1/p(1)
A=step4(10)+(ddt(10)./Qest);
QMMFS2=-ddt./(step4-A);

Step5=step3-step4;
Intri= -ddt./(Step5-A);

%Absolute attenuation values by a linear fitting of 8 receivers. 1<=i<=8

% x=tp(1:dp1);
% y=step4(1:dp1);
% p = polyfit(x,y,1);
% 
% x2=tp(1:dp2);
% y2=step4(1:dp2);
% p2 = polyfit(x2,y2,1);
% 
% x3=tp(1:dp3);
% y3=step4(1:dp3);
% p3 = polyfit(x3,y3,1);
% 
% x4=tp(1:dp4);
% y4=step4(1:dp4);
% p4 = polyfit(x4,y4,1);
% 
% Qest1=-1/p(1); 
% Qest2=-1/p2(1); 
% Qest3=-1/p3(1); 
% Qest4=-1/p4(1); 
% 
% Af1=step4(dp1)+(1/Qest1)*(tp(dp1));QMMFS1=(-tp./(step4-Af1));
% Af2=step4(dp2)+(1/Qest2)*(tp(dp2));QMMFS2=(-tp./(step4-Af2));
% Af3=step4(dp3)+(1/Qest3)*(tp(dp3));QMMFS3=(-tp./(step4-Af3));
% Af4=step4(dp4)+(1/Qest4)*(tp(dp4));QMMFS4=(-tp./(step4-Af4));
% 
% QMMFSF=(mean([QMMFS1' QMMFS2' QMMFS3' QMMFS4']'))';
% 
% figure;
% plot(QMMFS3(1:length(QMMFS3)),z(1:length(QMMFS3)));flipy; xlabel('Q'); ylabel('depth (m)'); title('QT estimation');
