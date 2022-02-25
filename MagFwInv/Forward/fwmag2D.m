% Forward model for Magnetometry in 2D
% It uses a grid with mag susc contrasts as input to calculate the magnetic
% profile on top of the model in nT. 
%
% Inputs:
% model: mag. Susceptibility contrast model of dimensions nx*nz
% H: Total magnetic field magnitude
% st: Geographical azimuth of the profile (Liu, 2015)
% dec: Declination of the magnetic field in degrees in the area
% I: Inclination of the magnetic field in degrees in the area
% nobs: number of data points to calculate
% dx: distance step along x
% dz: distance step along the vertical direction
% xaxis: axis of distance along x
% zaxis: axis of distance along the vertical direction
% flag: 1 for plotting
% 
% Outputs:
% d: magnetic profile at 'xobs' positions in nT
% xobs: positions of magnetic points along profile
% G: transfer function for the forward problem
% 
% J. B. Tary, Universidad de los Andes, June 18, 2018
% C. V. Barrera López, Universidad de los Andes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d,xobs,G] = fwmag2D(model,H,st,dec,nobs,I,dx,dz,xaxis,zaxis,flag)

nx = size(model,2); nz = size(model,1);
maxdist = max(xaxis)+(0.5*dx);
dobs = maxdist/(nobs-1); % dobs has 'nobs-1' points spanning 'maxdist' distance
xobs = 0:dobs:maxdist; % horizontal distance in x direction for obs

% Create a matrix G with the transfer function for this problem
G = zeros(nobs,(nz)*(nx));
for i=1:nobs
    xx = xaxis-xobs(i); % Distance relative to the observation points
    [X,Z] = meshgrid(xx,zaxis); % Grid of the x,z coordinates
    Y = 0; % Extra dimension for 3D prism (not used for now)

    %Vertices of the prism (following Liu et al., 2013)
    x1=X+dx;
    x2=x1;
    x3=X-dx;
    x4=x3;
    z1=Z-dz;
    z2=Z+dz;
    %Following Liu (2013)
    E=log((((z2.*z2)+(x2.*x2)).*((z1.*z1)+(x3.*x3)))./...
        (((z1.*z1)+(x1.*x1)).*((z2.*z2)+(x4.*x4))));
    Eij =reshape(E',[1,(nx)*(nz)]);
    F=atan((2.*dx.*z1)./((z1.*z1)+(x1-dx).^2-(dx.*dx)))-...
        atan((2.*dx.*z2)./((z2.*z2)+(x2-dx).^2-(dx.*dx)));
    Fij =reshape(F',[1,(nx)*(nz)]);

    % Calculate the 'G' matrix following formula of Magnetic Field
    str = st*pi/180;
    Decr = dec*pi/180;
    A = str-Decr; %Angle between profile and Geomagnetic north
    Ir = -I*pi/180;
    DH=2*H*((sin(Ir)*(Eij/2))-cos(Ir)*Fij);
    DZ=2*H*((cos(Ir)*(Eij/2))+sin(Ir)*Fij); 
    
    G(i,:) = DH*cos(Ir)*cos(A)+DZ*sin(Ir);
    
    clear xx X Z rr rij* zij
end

G(isnan(G)==1) = 0; % For rij=0, G goes to infinity and Matlab returns NaN
% Calculate the predicted data using your model and 'G'
mij = reshape(model',[1,(nx)*(nz)]);
d = G*mij';
if flag == 1
    figure;
    subplot(212); imagesc(xaxis,zaxis,model);
    xlabel('Distance along x (m)'); ylabel('Depth (m)')
    subplot(211); plot(xobs,d,'k');grid
    xlabel('Distance along x (m)'); ylabel('Magnetic anomaly (nT)')
end
