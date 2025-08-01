function [ug,vg] = clc_geocurrent(x,y,f,zeta,rho,z,type);
    %%%%%%%%%%%%%
    % x,y 坐标
    % zeta,自由海面高度
    % rho,密度
    % z 一维的深度坐标
    % type 选择从海表还是海底积分
    % 'top'从海表积分 直接使用地转流进行计算是没用Boussinesq近似的
    % 'bottom'从最深处积分 热成风向上积分默认Boussinesq近似，公式推导中带的
    %%%%%%%%%%%%%
    rho_r=1025;g=9.8;
    if strcmp(type,'top')
        ug=zeros(size(rho));
        vg=zeros(size(rho));
        p=zeros(size(rho));
        p(:,:,end)=rho(:,:,end).*g.*zeta;
        [dpdx,dpdy]=model_gradient(x,y,p(:,:,end));
        dh=z(2:end)-z(1:end-1);
        ug(:,:,end)=-dpdy./rho(:,:,end)./f;
        vg(:,:,end)=dpdx./rho(:,:,end)./f;
        for kk=length(z)-1:-1:1
            p(:,:,kk)=(rho(:,:,kk+1)+rho(:,:,kk))./2.*g.*dh(kk)+p(:,:,kk+1);
            [dpdx,dpdy]=model_gradient(x,y,p(:,:,kk));
            ug(:,:,kk)=-dpdy./rho(:,:,kk)./f;
            vg(:,:,kk)=dpdx./rho(:,:,kk)./f;
        end
    elseif strcmp(type,'bottom')
        ug=zeros(size(rho));
        vg=zeros(size(rho));
        dh=z(2:end)-z(1:end-1);
        for kk=2:length(z)
            [drhodx,drhody]=model_gradient(x,y,(rho(:,:,kk)+rho(:,:,kk-1))./2);
            drhodx(isnan(drhodx))=0;drhody(isnan(drhody))=0;
            ug(:,:,kk)=g.*drhody./f./rho_r.*dh(kk-1)+ug(:,:,kk-1);
            vg(:,:,kk)=-g.*drhodx./f./rho_r.*dh(kk-1)+vg(:,:,kk-1);
        end
    else
        error(message('type字符出错'));
    end
    
end