function [zeta,temp_output,salt_output,u_output,v_output,w_output,akv_output,akt_output,visc3d_output] = get_data_read_interp(grdname,filelist,filenum,N,theta_s,theta_b,hc,vtransform,z,h_deepest)
    %%%%%%%%%%%%%
    % gradname 网格文件路径
    % filelist 批量读取文件列表，主程序中获取方法，例：
    %       rpath='I:\ROMS_WP22_SCS_zheng_0.5km\avg\';
    %       list=dir(fullfile(rpath,'*avg*.nc.2'));
    % filenu 获取文件列表中第几个到第几个的数据
    % N，theta_s ，theta_b ，hc，vtransform S坐标参数
    % z 需要插值到的深度
    % h_deepest 大于多少的深度不予插值，这个是节省计算资源用的
    %%%%%%%%%%%%%
    
    h=ncread(grdname,'h',[598 194],[501 536]);f=ncread(grdname,'f',[598 194],[501 536]);
    x_rho= ncread(grdname,'x_rho',[598 194],[501 536]);   y_rho= ncread(grdname,'y_rho',[598 194],[501 536]);   % X/Y rho  :km
%     h=ncread(grdname,'h');f=ncread(grdname,'f');
%     x_rho= ncread(grdname,'x_rho');   y_rho= ncread(grdname,'y_rho');   % X/Y rho  :km
    step1=0;
    for ii=filenum
        step1=step1+1;
        disp(['stepnum: ',num2str(step1)]);
        disp(datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        filename=[filelist(ii).folder,'\',filelist(ii).name];
        temp=ncread(filename,'temp',[598 194 1 1],[501 536 inf 1]);
        salt=ncread(filename,'salt',[598 194 1 1],[501 536 inf 1]);
        u_u=ncread(filename,'u',[597 194 1 1],[502 536 inf 1]);
        v_v=ncread(filename,'v',[598 193 1 1],[501 537 inf 1]);
        for iii=1:size(u_u,3)
            u(:,:,iii)=v2rho_2d(u_u(:,:,iii));
            v(:,:,iii)=u2rho_2d(v_v(:,:,iii));
        end
        w=ncread(filename,'w',[598 194 1 1],[501 536 inf 1]);
        akv1=ncread(filename,'AKv',[598 194 1 1],[501 536 inf 1]);
        akv=(akv1(:,:,1:end-1)+akv1(:,:,2:end))./2;
        akt1=ncread(filename,'AKt',[598 194 1 1],[501 536 inf 1]);
        akt=(akt1(:,:,1:end-1)+akt1(:,:,2:end))./2;
        visc3d=ncread(filename,'visc3d',[598 194 1 1],[501 536 inf 1]);
        zeta=ncread(filename,'zeta',[598 194 1],[501 536 1]);
        %%%%
%         temp=ncread(filename,'temp');
%         salt=ncread(filename,'salt');
%         u_u=ncread(filename,'u');
%         v_v=ncread(filename,'v');
%         for iii=1:size(u_u,3)
%             u(:,:,iii)=v2rho_2d(u_u(:,:,iii));
%             v(:,:,iii)=u2rho_2d(v_v(:,:,iii));
%         end
%         w=ncread(filename,'w');
%         akv1=ncread(filename,'AKv');
%         akv=(akv1(:,:,1:end-1)+akv1(:,:,2:end))./2;
%         akt1=ncread(filename,'AKt');
%         akt=(akt1(:,:,1:end-1)+akt1(:,:,2:end))./2;
%         visc3d=ncread(filename,'visc3d');
%         zeta=ncread(filename,'zeta');
        
        % 计算z_rho
        z_rho=zlevs(h,zeta,theta_s,theta_b,hc,N,'r',vtransform);

        % 将theta插值到z
        parfor k = 1:size(x_rho,2)
            tempSlice = zeros(size(y_rho,1), 1, length(z));
            saltSlice = zeros(size(y_rho,1), 1, length(z));
            uSlice = zeros(size(y_rho,1), 1, length(z));
            vSlice = zeros(size(y_rho,1), 1, length(z));
            wSlice = zeros(size(y_rho,1), 1, length(z));
            akvSlice = zeros(size(y_rho,1), 1, length(z));
            aktSlice = zeros(size(y_rho,1), 1, length(z));
            visc3dSlice = zeros(size(y_rho,1), 1, length(z));

            for m = 1:size(y_rho,1)
                z_rho_use=squeeze(z_rho(:,m,k));z_rho_use=z_rho_use(z_rho_use>h_deepest);z_num=length(z_rho_use);
                tempSlice(m,:) = interp1(z_rho_use, squeeze(temp(m,k,end-z_num+1:end)), z, 'linear','extrap');
                saltSlice(m,:) = interp1(z_rho_use, squeeze(salt(m,k,end-z_num+1:end)), z, 'linear','extrap');
                uSlice(m,:) = interp1(z_rho_use, squeeze(u(m,k,end-z_num+1:end)), z, 'linear','extrap');
                vSlice(m,:) = interp1(z_rho_use, squeeze(v(m,k,end-z_num+1:end)), z, 'linear','extrap');
                wSlice(m,:) = interp1(z_rho_use, squeeze(w(m,k,end-z_num+1:end)), z, 'linear','extrap');
                akvSlice(m,:) = interp1(z_rho_use, squeeze(akv(m,k,end-z_num+1:end)), z, 'linear','extrap');
                aktSlice(m,:) = interp1(z_rho_use, squeeze(akt(m,k,end-z_num+1:end)), z, 'linear','extrap');
                visc3dSlice(m,:) = interp1(z_rho_use, squeeze(visc3d(m,k,end-z_num+1:end)), z, 'linear','extrap');
            end

            temp_output(:,k,:,step1) = tempSlice;
            salt_output(:,k,:,step1) = saltSlice;
            u_output(:,k,:,step1) = uSlice;
            v_output(:,k,:,step1) = vSlice;
            w_output(:,k,:,step1) = wSlice;
            akv_output(:,k,:,step1) = akvSlice;
            akt_output(:,k,:,step1) = aktSlice;
            visc3d_output(:,k,:,step1) = visc3dSlice;
        end   
    end
    % mask 地形
    for m = 1:size(x_rho,1)
        parfor k = 1:size(x_rho,2)
            mask_z(m,k,:)=z<-h(m,k);
        end
    end
    mask_z=repmat(mask_z,[1 1 1 step1]);
    temp_output(mask_z)=nan;
    salt_output(mask_z)=nan;
    u_output(mask_z)=nan;
    v_output(mask_z)=nan;
    w_output(mask_z)=nan;
    akv_output(mask_z)=nan;
    akt_output(mask_z)=nan;
    visc3d_output(mask_z)=nan;
    
end
