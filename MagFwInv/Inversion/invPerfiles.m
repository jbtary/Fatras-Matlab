%% Leer Perfiles de Anomalias
clear all; 
clc

P1 = xlsread('P1.xlsx');
P2 = xlsread('P2.xlsx');
P3 = xlsread('P3.xlsx');
P4 = xlsread('P4.xlsx');


%% Distancias y Magnitudes
%Tiempos
tP1 = P1(:,5);tP2 = P2(:,5);tP3 = P3(:,5);tP4 = P4(:,5);
%Distancias
xP1 = P1(:,7);xP2 = P2(:,7);xP3 = P3(:,7);xP4 = P4(:,7);
%Alturas
zP1 = P1(:,4);zP2 = P2(:,4);zP3 = P3(:,4);zP4 = P4(:,4);
%Anomalias
mgP1 = P1(:,6);mgP2 = P2(:,6);mgP3 = P3(:,6);mgP4 = P4(:,6);

figure
subplot(2,2,1);plot(xP1,mgP1,'g*-'); xlabel('Distance (m)');ylabel('Magnetic Anomaly (nT)')
subplot(2,2,2);plot(xP2,mgP2,'r*-'); xlabel('Distance (m)');ylabel('Magnetic Anomaly (nT)')
subplot(2,2,3);plot(xP3,mgP3,'k*-'); xlabel('Distance (m)');ylabel('Magnetic Anomaly (nT)')
subplot(2,2,4);plot(xP4,mgP4,'b*-'); xlabel('Distance (m)');ylabel('Magnetic Anomaly (nT)')

%Time
% figure (2)
% subplot(2,2,1);plot(tP1,mgP1,'g*-')
% subplot(2,2,2);plot(tP2,mgP2,'r*-')
% subplot(2,2,3);plot(tP3,mgP3,'k*-')
% subplot(2,2,4);plot(tP4,mgP4,'b*-')


%% least-squares solution
d= mgP1;
xProf = xP1;
nx=81;
nz=41;
I=29.2142;
H=31076.2; %Tot Mag field (26/07/2018) in the prospecting area.
dec = -7.4713; %Geomagnetic declination
st=0;
maxdept=60;
[G] = constG(d,xProf,I,H,st,dec,nx,nz,maxdept);

mest_ij = inv(G'*G)*transpose(G)*d;
dest = G*mest_ij;

mest = reshape(mest_ij,[(nx),(nz)]); 
mest = mest';

e = d-dest;
Erms = (e'*e)/length(d);
figure; 
subplot(211); hold on; title('LSQ method')
plot(xobs,d); plot(xobs,dest,'r--'); grid;
ylabel('Magnetic anomaly (nT)'); legend('Observed Data', 'Predicted data')
 
subplot(212); imagesc(xaxis,zaxis,mest);
xlabel('Distance along x (m)'); ylabel('Depth (m)');c = colorbar('southoutside');c.Label.String = 'Magnetic Susceptibility(SI)';


%% Damped Least-squares solution
% 10 - Use the damped least-squares solution to find an estimation of the 
% model from the G and d arrays (from Sacchi, M. D.,2006) 
clear mest* dest* mu Erms modnorm

mu_min=log10(1e6);mu_max = log10(1e9);Nmu=11;
M = length(G); 
% N = length(d);
% I=eye(N);
for k=1:Nmu
   mu(k)=mu_min+(mu_max-mu_min)*(k-1)/(Nmu-1);
   mu(k) = 10^mu(k);
%    m_est = G'*((G*G'+mu(k)*I)\d);
   m_est = inv(G'*G + mu(k)*eye(M))*G'*d;   
   dp = G*m_est;
   misfit(k)=(dp-d)'*(dp-d);
   modnorm(k) = m_est'*m_est;
   
   disp(['Loop #' num2str(k)])

end

%Fancy plot
figure(1);clf;
plot(modnorm,misfit,'*');hold on
plot(modnorm,misfit);

for k=Nmu:-1:8
    say=strcat('    \mu=',num2str(mu(k)));
    text(modnorm(k),misfit(k),say);
end
xlabel('Model Norm');ylabel('Misfit');grid
title('L-Curve')

%%
muSel = mu(8);
mest_ij = G'*((G*G'+muSel*I)\d);
% mest_ij =inv(G'*G + mu(3)*eye(M))*G'*d;
mest = reshape(mest_ij,[(nx),(nz)]); mest = mest';

% Calculate the predicted data using your estimated model and 'G'
dest = G*mest_ij;
e = d-dest;
Erms = (e'*e)/length(d);

figure; 
subplot(211); hold on; title('D-LSQ method')
plot(xobs,d); plot(xobs,dest,'r--');grid;
ylabel('Magnetic anomaly (nT)'); legend('Observed Data', 'Predicted data')
 
subplot(212); imagesc(xaxis,zaxis,mest);
xlabel('Distance along x (m)'); ylabel('Depth (m)');c = colorbar('southoutside');c.Label.String = 'Magnetic Susceptibility(SI)';


%% LSQ Weighted Minimum Norm (Following Sacchi, M. D. (2006) )
W = convmtx([1,-1],M);
W1 = W(1:M,1:M);
m_est = (G'*G+muSel*W1'*W1)\(G'*d);
dest=G*m_est;

mest = reshape(m_est,[nx,nz]); mest = mest';

figure; 
subplot(211); hold on; title('LSQWMN method')
plot(xobs,d); plot(xobs,dest,'r--'); grid;
ylabel('Magnetic anomaly (nT)'); legend('Observed Data', 'Predicted data')
 
subplot(212); imagesc(xaxis,zaxis,mest);
xlabel('Distance along x (m)'); ylabel('Depth (m)');c = colorbar('southoutside');c.Label.String = 'Magnetic Susceptibility(SI)';

%Error
modnorm = (dest'.*dest);
e = G*m_est-d;
ERMS = sqrt((e'*e)/length(d));

%% Conjugate gradient solution to the inversion problem
m0 = ones(length(G),1);

[mest_ij,diag] = cgls(G,m0,d,200,1e-20);
mest = reshape(mest_ij,[(nx),(nz)]); mest = mest';

% Calculate the predicted data using your estimated model and 'G'
dest = G*mest_ij;

figure; 
subplot(211); hold on; title('Conjugate Gradient solution')
plot(xobs,d); plot(xobs,dest,'r--'); grid;
ylabel('Magnetic anomaly (nT)'); legend('Observed Data', 'Predicted data')
 
subplot(212); imagesc(xaxis,zaxis,mest);
xlabel('Distance along x (m)'); ylabel('Depth (m)');c = colorbar('southoutside');c.Label.String = 'Magnetic Susceptibility(SI)';

e = G*mest_ij-d;
ERMS = sqrt((e'*e)/length(d));

