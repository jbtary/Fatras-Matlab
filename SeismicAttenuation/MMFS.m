function [QMMFS,QMMFS1,step4,A,step3,bet,step2,step1]=MMFS(yy,ddt,x,dep,realQ,testQ,z,Q)%,dp1,dp2,dp3,dp4)

bet=yy./(pi.*x');

step1=mean(bet);

step2=(median((bet-step1)'))';

step3=bet-step2;

step4=median(step3);

A=step4(dep)+(ddt(dep)./realQ);
A1=step4(dep)+(ddt(dep)./testQ);
%A2=step4(dep)+(ddt(dep)./testQ2);
QMMFS=-ddt./(step4-A);
QMMFS1=-ddt./(step4-A1);
%QMMFS2=-ddt./(step4-A2);


figure;
plot((1./QMMFS(1:(length(QMMFS)))),z(1:(length(QMMFS))),'o');flipy; xlabel('1/Q'); ylabel('depth (m)'); title('DCI+QMMFS estimation');
hold on;plot(1./QMMFS1(1:length(QMMFS1)),z(1:length(QMMFS1)),'o');
hold on;plot((1./Q(1:(length(QMMFS1)))),z(1:(length(QMMFS1))));
%hold on;plot((1./QMMFS2(1:(length(QMMFS2)))),z(1:(length(QMMFS2))),'o');
legend('Absulute Q by MMFS.','Real Q.','Relative Q by MMFS (Q = 30).')

Sca=step3-step4;


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
