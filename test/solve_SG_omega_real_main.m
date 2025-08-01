clear all; warning off all;  close all;  warning('off'); warning off;
cd E:\SG-Omega
%read grid
grdname='I:\ROMS_WP22_SCS_zheng_1.5km_fromrst480\preprocessing flie\roms_grd.nc.2';
% 截取区域，不一定全要啊
h=ncread(grdname,'h',[598 194],[501 536]);f=ncread(grdname,'f',[598 194],[501 536]);
lon_rho= ncread(grdname,'lon_rho',[598 194],[501 536]);   lat_rho= ncread(grdname,'lat_rho',[598 194],[501 536]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[598 194],[501 536]);   y_rho= ncread(grdname,'y_rho',[598 194],[501 536]);   % X/Y rho  :km

N= 60; theta_s= 7; theta_b= 2; hc= 100;vtransform= 2.;
rho_r=1025;g=9.8;

z=[-300:5:0];
h_deepest=-400;
rho_r=1025;g=9.8;

rpath1='I:\ROMS_WP22_SCS_zheng_0.5km\avg\';
filelist=dir(fullfile(rpath1,'*avg*.nc.2'));
filenum=24:26;
[zeta1,temp,salt,u,v,w,akv,akt,visc3d]=get_data_read_interp(grdname,filelist,filenum,N,theta_s,theta_b,hc,vtransform,z,h_deepest);

rho=sw_dens0(salt,temp);
for ii=1:size(rho,4);
    [ug(:,:,:,ii),vg(:,:,:,ii)] = clc_geocurrent(x_rho,y_rho,f,zeta1,rho(:,:,:,ii),z,'bottom');
end
dt=2;
[Q]=clc_Qvector(x_rho,y_rho,z,dt,f,rho,u,v,ug,vg,akv,akt,visc3d);

w_omega=solve_SG_omega(x_rho,y_rho,z,rho(:,:,:,2:end-1),Q,f,1.5,100,1e-20);

% [dudx,dudy]=model_gradient(x_rho,y_rho,u(:,:,:,2));
% [dvdx,dvdy]=model_gradient(x_rho,y_rho,v(:,:,:,2));
% zeta=(dvdx-dudy)./f;
figure
pcolor(lon_rho,lat_rho,w_omega(:,:,end-3));shading interp
caxis([-5e-4 5e-4])
figure
pcolor(lon_rho,lat_rho,w(:,:,end-3));shading interp
caxis([-5e-4 5e-4])
