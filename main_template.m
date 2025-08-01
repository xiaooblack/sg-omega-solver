% main
% 这份代码是为怎样正确使用工具包中函数进行计算而写的主程序模板
clear all; warning off all;  close all
%% step 1：获取数据
% 我们常用的三维数据或者模式数据，需要将其转化为单位为m的网格
% 这里示例是roms的模式数据
% utils中的 get_data_read_interp 函数专门为处理模式数据而设计
grdname='I:\ROMS_WP22_SCS_zheng_1.5km_fromrst480\preprocessing flie\roms_grd.nc.2';
h=ncread(grdname,'h',[598 194],[501 536]);f=ncread(grdname,'f',[598 194],[501 536]);
lon_rho= ncread(grdname,'lon_rho',[598 194],[501 536]);   lat_rho= ncread(grdname,'lat_rho',[598 194],[501 536]);   % X/Y rho  :km
x_rho= ncread(grdname,'x_rho',[598 194],[501 536]);   y_rho= ncread(grdname,'y_rho',[598 194],[501 536]);   % X/Y rho  :km
N= 60; theta_s= 7; theta_b= 2; hc= 100;vtransform= 2.; %
rho_r=1025;g=9.8;
z=[-300:5:0];% 垂向网格自定义，可以不均匀
h_deepest=-400;% 选取插值所用的最深深度
rho_r=1025;g=9.8;
rpath1='I:\ROMS_WP22_SCS_zheng_0.5km\avg\';
filelist=dir(fullfile(rpath1,'*avg*.nc.2'));% 模式数据存放文件夹
filenum=24:26;% nc文件夹中的第n到m个文件，三个以上才能进行计算，因为计算非地转流局地时间变化的Q-vector采用的是中心插值
[zeta1,temp,salt,u,v,w,akv,akt,visc3d]=get_data_read_interp(grdname,filelist,filenum,N,theta_s,theta_b,hc,vtransform,z,h_deepest);
%% step 2: 计算地转流
% utils中的 clc_geocurrent 函数为计算地转流而设计，配备了两种计算方法，从底层用热成风积分以及从海表计算正压与斜压分量
% zeta1 变量为海表高度，如果采用底层热成风积分到海表的方法，这个zeta1给个0数组即可
% f为2d科氏力参数
rho=sw_dens0(salt,temp); % 获取密度
for ii=1:size(rho,4);
    [ug(:,:,:,ii),vg(:,:,:,ii)] = clc_geocurrent(x_rho,y_rho,f,zeta1,rho(:,:,:,ii),z,'bottom');
end
%% step 3：计算Q vector
dt=2;% 给予时间层之间的时间间隔，单位为小时
[Q]=clc_Qvector(x_rho,y_rho,z,dt,f,rho,u,v,ug,vg,akv,akt,visc3d);
% 输出的Q为结构数组，且对应第n+1至m-1个文件，两端不是中心插值就不要了
%% step 4：离散方程并使用超松弛迭代法进行计算
w_omega=solve_SG_omega(x_rho,y_rho,z,rho(:,:,:,2:end-1),Q,f,1.5,100,1e-20);
% 注意这里的rho对应为第n+1至m-1层
% 可以手动修改solve_SG_omega函数中的Q vector，用于判断不同强迫的贡献
