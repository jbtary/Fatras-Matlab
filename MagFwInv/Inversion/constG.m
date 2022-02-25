%Construct the mathematical description of the Magnetic 2D problem
%following Liu (2013).
%
%Inputs:
% d= Measurements
% xProf = Positions of the measurements (in m)
% I: inclination of the field in degrees
% H: Total magnetic field magnitude
% st: Geographical azimuth of the profile (Liu, 2015)
% dec: Geomagnetic declination in the area
% nx: number of prisms along x
% nz: number of prisms along the vertical direction
% maxdept: Max dept of the model
%
%Outputs:
% G: mathematical description of the model
%
%Carol Vanessa Barrera López, Universidad de los Andes, Nov 14, 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G] = constG(d,xProf,I,H,st,dec,nx,nz,maxdept)
nobs=length(d);
maxdist=round(max(xProf));
dx = maxdist/(nx-1); dz = maxdept/(nz-1);
xaxis = 0+dx/2:dx:maxdist;
zaxis = 0+dz/2:dz:maxdept;
dobs = maxdist/(nobs-1); % dobs has 'nobs-1' points spanning 'maxdist' distance
xobs = 0:dobs:maxdist; % horizontal distance in x direction for obs

% Update the number of nodes
nx = nx - 1; nz = nz - 1;
% Create a matrix G with the transfer function for this problem
G = zeros(nobs,(nz)*(nx));
for i=1:nobs
    xx = xaxis-xobs(i); % Distance relative to the observation points
    [X,Z] = meshgrid(xx,zaxis); % Grid of the x,y coordinates
    Y = 0; % Extra dimension for 3D prism (not used for now)

    x1=X+dx;
    x2=x1;
    x3=X-dx;
    x4=x3;
    z1=Z-dz;
    z2=Z+dz;
    E=log((((z2.*z2)+(x2.*x2)).*((z1.*z1)+(x3.*x3)))./...
        (((z1.*z1)+(x1.*x1)).*((z2.*z2)+(x4.*x4))));
    Eij =reshape(E',[1,(nx)*(nz)]);
    F=atan((2.*dx.*z1)./((z1.*z1)+(x1-dx).^2-(dx.*dx)))-...
        atan((2.*dx.*z2)./((z2.*z2)+(x2-dx).^2-(dx.*dx)));
    Fij =reshape(F',[1,(nx)*(nz)]);

    % Calculate the 'G' matrix following formula of vertical rectangular
    % prism
    str = st*pi/180;
    Decr = dec*pi/180;
    A = str-Decr; %Angle between profile and Geomagnetic north
    Ir = -I*pi/180;
    DH=2*H*((sin(Ir)*(Eij/2))-cos(Ir)*Fij);
    DZ=2*H*((cos(Ir)*(Eij/2))+sin(Ir)*Fij); 
    
    G(i,:) = (DH*cos(Ir)*cos(A)+DZ*sin(Ir));
    
    clear xx X Z rr rij* zij DH DZ E Eij F Fij st str dec Decr A Ir H x1 x2 x3 x4
end

G(isnan(G)==1) = 0; % For rij=0, G goes to infinity and Matlab returns NaN
end
