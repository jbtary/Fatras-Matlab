% Build customizable density model for gravity data modeling

function [model,mij,xaxis,zaxis,dx,dz] = modbuildgrav(nx,nz,maxdist,maxdept)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   nx: number of grid nodes along the x direction
%   nz: number of grid nodes along the vertical direction
%   maxdist: maximum distance along the x direction (thought to be in
%   meters but don't have to be). Distance along x will then be from 0 to
%   maxdist.
%   maxdept: maximum depth of the model (thought to be in meters as well
%   but don't have to be). Distance along the vertical direction will then
%   be from 0 to maxdept.
%   
% Outputs:
%   model: final matrix (nz-1)*(nx-1) of the model
%   mij: final model in vector form 1*((nx-1)*(nz-1)) for inversion
%   xaxis: axis of distance along x
%   xaxis: axis of distance along the vertical direction
%   dx: distance step along x
%   dz: distance step along the vertical direction
% 
% Example:
% [m,mij,xaxis,zaxis,dx,dz] = modbuildgrav(101,51,200,50);
% 
% J.B. Tary, Universidad de los Andes, 18/04/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a 2D grid of nodes with a background density
dx = maxdist/(nx-1); dz = maxdept/(nz-1);
xaxis = 0+dx/2:dx:maxdist;
zaxis = 0+dz/2:dz:maxdept;
% Update the number of nodes
nx = nx - 1; nz = nz - 1;
model = zeros(nz,nx);

[X,Z] = meshgrid(xaxis,zaxis); % Grid of the x,y coordinates
zij = reshape(Z',[1,nx*nz]); % Reshape depth 'Z' to a column vector
xij = reshape(X',[1,nx*nz]); % Reshape depth 'Z' to a column vector

% In the model, put an/some area/s with different density anomalies
top = 1;
while top ~= 0
    figure(87);
    imagesc(xaxis,zaxis,model); 
    xlabel('Distance along x'); ylabel('Depth')
    set(gca,'XTick',0:dx:maxdist); set(gca,'XTickLabel',[])
    set(gca,'YTick',0:dz:maxdept); set(gca,'YTickLabel',[])
    colorbar; grid on
    
    disp('Create a polygon corresponding to a zone of different density')
    h = impoly;
    pos = getPosition(h);
    [in,on] = inpolygon(xij,zij,pos(:,1),pos(:,2));
    inon = in | on;
    
    prompt = 'What is the value of density contrast you want inside this polygon (ex: 500)? ';
    density = input(prompt);
    
    idx = find(inon(:)); % Linear indices of 'inon' points
    
    mij = reshape(model',[1,nx*nz]);
    mij(idx) = density;
    mij = reshape(mij,[nx,nz]); mij = mij';
    model = mij;
    
    figure(87);
    imagesc(xaxis,zaxis,model); 
    xlabel('Distance along x'); ylabel('Depth')
    set(gca,'XTick',0:dx:maxdist); set(gca,'XTickLabel',[])
    set(gca,'YTick',0:dz:maxdept); set(gca,'YTickLabel',[])
    colorbar; grid on
    
    prompt = 'Press 0 if you are satisfied with your model, other numbers and you continue. ';
    top = input(prompt);
    
    clear density prompt h pos in on inon mij idx
end
close figure 87

mij = reshape(model',[1,nx*nz]); % Model as a vector to be used in the inversion
% Plotting
imagesc(xaxis,zaxis,model); colorbar
xlabel('Distance along x (m)'); ylabel('Depth (m)'); title('Final model')
