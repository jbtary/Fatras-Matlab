%% Building the Forward Problem
clear 
clc
% Make a 2D grid of nodes with a background density
% In the model, put an/some area/s with different density anomalies
nx=81;
nz=51;
maxdist=400;
maxdept=100;
[model,mij,xaxis,zaxis,dx,dz] = modbuildgrav(nx,nz,maxdist,maxdept);

figure;
subplot(212); imagesc(xaxis,zaxis,model); 
xlabel('Distance along x (m)'); ylabel('Depth (m)')

%% Prueba (Liu, 2013)
nobs=601;
I =29.2142; %(26/07/2018) in the prospecting area.
H=31076.2; %(26/07/2018) in the prospecting area.
dec = -7.4713; %Geomagnetic declination
% load Tequendama
% model= modCorte; %Model of the Tequendama mine
[d,xobs,G] = fwmag2D(model,H,0,dec,nobs,I,dx,dz,xaxis,zaxis,0);

 figure;
 subplot(212); imagesc(xaxis,zaxis,model);
 xlabel('Distance along x (m)'); ylabel('Depth (m)')
 subplot(211); plot(xobs,d,'k');grid
 xlabel('Distance along x (m)'); ylabel('Magnetic anomaly (nT)')

 
