function [] = GravFunc( root,ax_ymin,ax_xmax )

%%
% Computes g_z and dg_z/dz for  a polygon using Talwani
% line integral
% Need to check direction of loop definition (clockwise!)

% This version reads in density as well from file and can
% model multiple polygons
% Load density model
% Number of polygons = Npolygons
% ith polygon has density rho(i) and nvertex(i) corners
root;

fid =fopen([root,'_mod.txt'],'r');
background_density=fscanf(fid,'%f',1);
npolygons=fscanf(fid,'%i',1);
for ipoly=1:npolygons
   rho(ipoly) = fscanf(fid,'%f',1);
   drho(ipoly) = rho(ipoly)-background_density;
   nvertex(ipoly) = fscanf(fid,'%f',1);
   for iv=1:nvertex(ipoly)
     xz=fscanf(fid,'%f',2);
     xp(ipoly,iv)=xz(1)*1000;
     zp(ipoly,iv)=xz(2)*1000;
   end
end


% Read in gravity data
fid =fopen([root,'_X.txt'],'r');
alln=fscanf(fid, '%f');
fclose(fid);
ndata = length(alln)/3;
fid = fopen([root,'_X.txt'],'r');
for i =1:ndata
  data = fscanf(fid,'%f',3);
  x_data(i)  = data(1);
  g_data(i)  = data(2);
  g_error(i) = data(3);
end
% Find maximum values of offset and Bouguer data
  xmax = 1.1*max(x_data);
  xmin = 1.1*min(x_data);
  ndata = length(x_data);

  g_max=0.9*max(g_data);
  g_min=1.1*min(g_data);

  
% Vector with colour flags for polygons
col=['g','b','c','r','m','y'];
G = 6.67e-11;
% No. of integration points on each polygon side

npts = 60;
xobs = x_data;

nx =length(xobs);
zmin = 1.2*max(max(zp));
npoints =length(xp);


% start loop over gravity measurement points
for ix=1:nx   
% Loop over polygons   
  for ipoly=1:npolygons
  % Loop over each side of polygon
    for iv =1:nvertex(ipoly)-1   
      dz = (zp(ipoly,iv+1)-zp(ipoly,iv))/npts;
      dx = (xp(ipoly,iv+1)-xp(ipoly,iv))/npts;
      side(ipoly,iv) = 0.0;
  % Integrate over side with npts steps - (Talwani's method integral) 
      for j = 1:npts
        zint = zp(ipoly,iv)+j*dz;
        if abs(zint) < 0.000001
           zint = 0.000001;
        end
        xint = xp(ipoly,iv)+j*dx;
        theta = atan(abs(xint-xobs(ix))/zint);
        if xobs(ix)>= xint
          theta = -theta;
        end
        side(ipoly,iv) = side(ipoly,iv) + theta*dz;
      end
    end   
  
    % Sum over the different sides of the polygon. This is needed since polygons
    % can have different numbers of sides. 
    dg(ipoly,ix) = 0.;
    for iv=1:nvertex(ipoly)-1
      dg(ipoly,ix) =(dg(ipoly,ix)-1.0e5*2*G*drho(ipoly)*side(ipoly,iv));
    end   
  end 
end

% Sum over contributions from different polygons (dg) to get Bouguer anomaly (gb)
for i=1:ndata
   gb(i)=0.0;
   for ipoly=1:npolygons
      gb(i)=gb(i)+dg(ipoly,i);
   end
end   

% Plot density model
ax1=[0 ax_xmax ax_ymin 0];
figure('pos',[0 0 434 600]); hold on

% Assign colour of polygon depending on density
subplot(211);
for ipoly=1:npolygons
    drho1=abs(rho(ipoly)-3000);
    drho2=abs(rho(ipoly)-3100);
    drho3=abs(rho(ipoly)-3200);
    drho4=abs(rho(ipoly)-3300);
    drho5=abs(rho(ipoly)-3400);
    drho6=abs(rho(ipoly)-3500);
    drho7=abs(rho(ipoly)-3600);
    if drho1<drho2 && drho1<drho3 && drho1<drho4 && drho1<drho5 && drho1<drho6 && drho1<drho7 
        c=[0.1 0 0.5];
    elseif drho2<drho1 && drho2<drho3 && drho2<drho4 && drho2<drho5 && drho2<drho6 && drho2<drho7 
        c=[0.1 0.5 1];
    elseif drho3<drho1 && drho3<drho2 && drho3<drho4 && drho3<drho5 && drho3<drho6 && drho3<drho7
        c=[0 1 0.1];
    elseif drho4<drho1 && drho4<drho2 && drho4<drho3 && drho4<drho5 && drho4<drho6 && drho4<drho7
        c=[1 1 0];
    elseif drho5<drho1 && drho5<drho2 && drho5<drho3 && drho5<drho4 && drho5<drho6 && drho5<drho7
        c=[1 0.7 0];
    elseif drho6<drho1 && drho6<drho2 && drho6<drho3 && drho6<drho4 && drho6<drho5 && drho6<drho7
        c=[1 0.6 0];
    elseif drho7<drho1 && drho7<drho2 && drho7<drho3 && drho7<drho5 && drho7<drho6 && drho7<drho7
        c=[0.9 0.2 0.1];
    else 
        c=[0.9 0 0];
    end
    fill(xp(ipoly,:),-1*zp(ipoly,:),c); hold on
end   

axis(ax1);
xlabel('Distance along the profile (m)'); ylabel('Depth (m)'); title(['Earth Density Model'])
hold off;

% Plot of Bouguer anomaly 
ax2=[0 ax_xmax -100 100];
ax2=subplot(212)
plot(xobs,gb,'-','LineWidth',2, 'color', '[0 0 0.3]');
hold on
%plot(x_data,g_data,'*'); %Si en Model_X columna 2 hubiesen datos medidos en campo para comparar
axis(ax2);
xlabel('Distance along the profile (m)'); ylabel('Bougeur Anomaly (miligals)'); title('Calculated Gravity Anomaly')
legend('Calculated Model', 'Location', 'southeast'); grid on
hold off


