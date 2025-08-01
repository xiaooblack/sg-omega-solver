function [Q]=clc_Qvector(x,y,z,dt,f,rho,u,v,ug,vg,akv,akt,visc3d);
    %%%%%%%%%%%%%
    % 计算Q vector 存成结构数组
    % ------输出---------
    % Q_tg 地转强迫Q vector
    % Q_tag 非地转强迫Q vector
    % Q_dag 地转和非地转有个交叉项，简化自广义Omega方程dag(因为这是半地转框架)
    % Q_th 海表动量输入，用动量混合系数算
    % Q_dm 海表热量输入，用温度混合系数算
    % Q_tr 多加的一项非地转流局地时间变化
    % ------输入---------
    % 对于输入的时间层只能计算 2到end-1的层
    % x,y 水平网格，二维
    % z 垂向坐标，一维
    % dt 时间层间隔,单位h
    % f 科氏力参数，二维
    % rho 密度
    % u,v 全流速
    % ug,vg 地转流
    % akv 垂向动量混合系数
    % akt 垂向温度混合系数
    % visc3d 水平动量混合系数
    %%%%%%%%%%%%%
    rho_r=1025;g=9.8;
    ua=u-ug;
    va=v-vg;
    b=-g.*rho./rho_r;
    step1=0;
    for ii=2:size(rho,4)-1
        step1=step1+1;
        [dugdx,dugdy]=model_gradient(x,y,ug(:,:,:,ii));
        [dvgdx,dvgdy]=model_gradient(x,y,vg(:,:,:,ii));
        [duadx,duady]=model_gradient(x,y,ua(:,:,:,ii));
        [dvadx,dvady]=model_gradient(x,y,va(:,:,:,ii));
        [dbdx,dbdy]=model_gradient(x,y,b(:,:,:,ii));
        duadz=model_gradient_z(ua(:,:,:,ii),z);
        dvadz=model_gradient_z(va(:,:,:,ii),z);
        dbdz=model_gradient_z(b(:,:,:,ii),z);
        [dudx,dudy]=model_gradient(x,y,u(:,:,:,ii));
        [dvdx,dvdy]=model_gradient(x,y,v(:,:,:,ii));
        [ddudxdx,~]=model_gradient(x,y,dudx);[~,ddudydy]=model_gradient(x,y,dudy);
        [ddvdxdx,~]=model_gradient(x,y,dvdx);[~,ddvdydy]=model_gradient(x,y,dvdy);
        dudz=model_gradient_z(u(:,:,:,ii),z);
        dvdz=model_gradient_z(v(:,:,:,ii),z);
        % Q_tg
        Q.Q_tgx(:,:,:,step1)=-2.*(dvgdx.*dbdy+dugdx.*dbdx);
        Q.Q_tgy(:,:,:,step1)=-2.*(dvgdy.*dbdy+dugdy.*dbdx);
        % Q_tag
        Q.Q_tagx(:,:,:,step1)=-(dvady.*dbdx+2.*duadx.*dbdx+dvadx.*dbdy);
        Q.Q_tagy(:,:,:,step1)=-(duady.*dbdx+2.*dvady.*dbdy+duadx.*dbdy);
        % Q_dag
        Q.Q_dagx(:,:,:,step1)=(dvadz.*dvgdy+duadz.*dvgdx).*repmat(f,[1 1 size(rho,3)]);
        Q.Q_dagy(:,:,:,step1)=(duadz.*dugdx+dvadz.*dugdy).*repmat(f,[1 1 size(rho,3)]).*-1;
        % Q_th
        Dvu=model_gradient_z(akv(:,:,:,ii).*dudz,z);
        Dhu=visc3d(:,:,:,ii).*(ddudxdx+ddudydy);
        Dvv=model_gradient_z(akv(:,:,:,ii).*dvdz,z);
        Dhv=visc3d(:,:,:,ii).*(ddvdxdx+ddvdydy);
        Q.Q_thx(:,:,:,step1)=model_gradient_z((Dvv+Dhv),z);
        Q.Q_thx(:,:,:,step1)=Q.Q_thx(:,:,:,step1).*repmat(f,[1 1 size(rho,3)]).*-1;
        Q.Q_thy(:,:,:,step1)=model_gradient_z((Dvu+Dhu),z);
        Q.Q_thy(:,:,:,step1)=Q.Q_thy(:,:,:,step1).*repmat(f,[1 1 size(rho,3)]);
        % Q_dm
        Dvb=model_gradient_z(akt(:,:,:,ii).*dbdz,z);
        [Q.Q_dmx(:,:,:,step1),Q.Q_dmy(:,:,:,step1)]=model_gradient(x,y,Dvb);
        % Q_tr
        duadz1=model_gradient_z(ua(:,:,:,ii-1),z);
        dvadz1=model_gradient_z(va(:,:,:,ii-1),z);
        duadz2=model_gradient_z(ua(:,:,:,ii+1),z);
        dvadz2=model_gradient_z(va(:,:,:,ii+1),z);
        Q.Q_trx(:,:,:,step1)=(dvadz2-dvadz1)./(dt*3600).*repmat(f,[1 1 size(rho,3)]);
        Q.Q_try(:,:,:,step1)=(duadz2-duadz1)./(dt*3600).*repmat(f,[1 1 size(rho,3)]).*-1;
    end
end