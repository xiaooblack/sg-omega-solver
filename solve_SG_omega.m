function w=solve_SG_omega(x,y,z,rho,Q,f, omega, max_iter, tolerance)
    %%%%%%%%%%%%%
    % 迭代法求解Omega方程,嫌弃收敛速度慢可以使用松弛迭代法
    % Omega方程是一个典型的斯特姆刘维尔问题，迭代法的前提是个PDE是椭圆型的，这样迭代才收敛
    % 很显然若是想满足是个椭圆型的方程，需要的是N2>=0;
    % ------输入---------
    % x,y 水平网格，二维
    % z 垂向坐标，一维
    % rho 密度，注意输入的时候把两端没算的时间层去掉
    % Q Q vector 结构数组
    % f 科氏力参数，二维数组
    %   omega 标量, 超松弛因子 (例如 1.8),1 就变成高斯迭代法
    %   max_iter 整数, 最大迭代次数 (例如 5000)
    %   tolerance 标量, 收敛阈值 (例如 1e-7)
    %%%%%%%%%%%%%
    rho_r=1025;g=9.8;
    for ii=1:size(rho,4)
        [dQ_tgxdx,~]=model_gradient(x,y,Q.Q_tgx(:,:,:,ii));
        [~,dQ_tgydy]=model_gradient(x,y,Q.Q_tgy(:,:,:,ii));
        [dQ_tagxdx,~]=model_gradient(x,y,Q.Q_tagx(:,:,:,ii));
        [~,dQ_tagydy]=model_gradient(x,y,Q.Q_tagy(:,:,:,ii));
        [dQ_dagxdx,~]=model_gradient(x,y,Q.Q_dagx(:,:,:,ii));
        [~,dQ_dagydy]=model_gradient(x,y,Q.Q_dagy(:,:,:,ii));
        [dQ_thxdx,~]=model_gradient(x,y,Q.Q_thx(:,:,:,ii));
        [~,dQ_thydy]=model_gradient(x,y,Q.Q_thy(:,:,:,ii));
        [dQ_dmxdx,~]=model_gradient(x,y,Q.Q_dmx(:,:,:,ii));
        [~,dQ_dmydy]=model_gradient(x,y,Q.Q_dmy(:,:,:,ii));
        [dQ_trxdx,~]=model_gradient(x,y,Q.Q_trx(:,:,:,ii));
        [~,dQ_trydy]=model_gradient(x,y,Q.Q_try(:,:,:,ii));
        % RHS 修改这里以调整不同的强迫
        RHS= dQ_trxdx+dQ_trydy + dQ_tgxdx+dQ_tgydy + dQ_tagxdx+dQ_tagydy + dQ_dagxdx+dQ_dagydy + dQ_thxdx+dQ_thydy + dQ_dmxdx+dQ_dmydy;
%        RHS= dQ_tgxdx+dQ_tgydy;
        % N2
        b=-g.*rho(:,:,:,ii)./rho_r;
        N2=model_gradient_z(b,z);
        N2(N2<0)=0;
        
        % 初始化和获取网格尺寸
        [nx, ny, nz] = size(N2);

        % 初始化垂向速度 w
        w1 = zeros(nx, ny, nz);
        % 边界条件 w=0 已经通过初始化满足 (Dirichlet条件)

        % 确保f是三维的，以简化后续计算
        if ismatrix(f)
            f = repmat(f, [1, 1, nz]);
        end
        f2 = f.^2;

        % 预计算网格间距
        % _f：i+1/2
        % _b: i-1/2
        
        % --- 垂直间距 (dz) ---
        % z 是一维的
        
        dz_f = zeros(nz); % 向前差分 z(k+1)-z(k)
        dz_b = zeros(nz); % 向后差分 z(k)-z(k-1)
        dz_f(2:nz-1) = z(2:end-1)-z(1:end-2);
        dz_b(2:nz-1)   = z(3:end)-z(2:end-1);
        % 边界上的间距可以设为与邻居相同，以避免在循环中出现零,反正边界点的也不计算用不到
        dz_f(nz) = dz_f(nz-1);dz_f(1) = dz_f(2);
        dz_b(nz) = dz_b(nz-1);dz_b(1) = dz_b(2);

        % --- 水平间距 (dx, dy) ---
        % x 和 y 是二维的
        dx_f = zeros(nx, ny); dy_f = zeros(nx, ny);
        dx_b = zeros(nx, ny); dy_b = zeros(nx, ny);
        % x-方向间距
        dx_f(2:nx-1, :) = x(2:end-1,:)-x(1:end-2,:);
        dx_b(2:nx-1, :)   = x(3:end,:)-x(2:end-1,:);
        dx_f(nx, :) = dx_f(nx-1, :);dx_f(1, :) = dx_f(2, :);
        dx_b(nx, :) = dx_b(nx-1, :);dx_b(1, :) = dx_b(2, :);

        % y-方向间距
        dy_f(:,2:ny-1) = y(:,2:end-1)-y(:,1:end-2);
        dy_b(:,2:ny-1) = y(:,3:end)-y(:,2:end-1);
        dy_f(:,ny) = dy_f(:,ny-1);dy_f(:,1) = dy_f(:,2);
        dy_b(:,ny) = dy_b(:,ny-1);dy_b(:,1) = dy_b(:,2);

        % 将水平间距扩展到三维
        dx_f = repmat(dx_f, [1, 1, nz]);
        dx_b = repmat(dx_b, [1, 1, nz]);
        dy_f = repmat(dy_f, [1, 1, nz]);
        dy_b = repmat(dy_b, [1, 1, nz]);

        % 预计算迭代系数
        % --- 邻居点系数 ---
        C_right=zeros(nx,ny,nz);C_left=zeros(nx,ny,nz);
        C_up=zeros(nx,ny,nz);C_down=zeros(nx,ny,nz);
        C_top=zeros(nx,ny,nz);C_bottom=zeros(nx,ny,nz);
        % 只计算内部点的系数即可
        % 内部点 (2:end-1)
        % i, j, k 索引是相对于整个数组的
        for k = 2:nz-1
            for j = 2:ny-1
                for i = 2:nx-1
                    % 垂直系数
                    dz_avg = 0.5 * (dz_f(k) + dz_b(k));
                    C_top(i,j,k) = (f2(i,j,k) / dz_avg) / dz_f(k);
                    C_bottom(i,j,k) = (f2(i,j,k) / dz_avg) / dz_b(k);

                    % 水平系数 (x-dir)
                    dx_avg = 0.5 * (dx_f(i,j,k) + dx_b(i,j,k));
                    N2_right = 0.5 * (N2(i+1,j,k) + N2(i,j,k));
                    N2_left = 0.5 * (N2(i,j,k)   + N2(i-1,j,k));
                    C_right(i,j,k) = (1 / dx_avg) * N2_right / dx_f(i,j,k);
                    C_left(i,j,k) = (1 / dx_avg) * N2_left / dx_b(i,j,k);

                    % 水平系数 (y-dir)
                    dy_avg = 0.5 * (dy_f(i,j,k) + dy_b(i,j,k));
                    N2_up = 0.5 * (N2(i,j+1,k) + N2(i,j,k));
                    N2_down = 0.5 * (N2(i,j,k)   + N2(i,j-1,k));
                    C_up(i,j,k) = (1 / dy_avg) * N2_up / dy_f(i,j,k);
                    C_down(i,j,k) = (1 / dy_avg) * N2_down / dy_b(i,j,k);
                end
            end
        end

        % --- 中心点系数的倒数 (用于除法) ---
        % 计算 1/C_center 来将除法变为乘法，效率更高。
        Inv_C_center = zeros(nx,ny,nz);
        C_center_val = (C_right + C_left + C_up + C_down + C_top + C_bottom);
        % 避免除以零，只在内部点计算
        Inv_C_center(2:nx-1, 2:ny-1, 2:nz-1) = 1 ./ C_center_val(2:nx-1, 2:ny-1, 2:nz-1);


        % 4. SOR 主循环
        fprintf('开始 SOR 迭代...\n');
        tic; % 开始计时
        for iter = 1:max_iter
            w_old = w1; % 保存上一步的结果以计算误差

            % 遍历所有内部网格点
            for k = 2:nz-1
                for j = 2:ny-1
                    for i = 2:nx-1
                        % 计算邻居项之和
                        sum_neighbors = ...
                            C_right(i,j,k) * w1(i+1, j, k)   + C_left(i,j,k) * w1(i-1, j, k) + ...
                            C_up(i,j,k) * w1(i, j+1, k)   + C_down(i,j,k) * w1(i, j-1, k) + ...
                            C_top(i,j,k) * w1(i, j, k+1)   + C_bottom(i,j,k) * w1(i, j, k-1);

                        % 计算雅可比估计值 w_jacobi
                        % w_jacobi = (RHS - sum_neighbors) / C_center
                        % 使用预计算的倒数: w_jacobi = (sum_neighbors - RHS) * Inv_C_center
                        w_jacobi = (sum_neighbors - RHS(i,j,k)) * Inv_C_center(i,j,k);

                        % 应用SOR更新公式
                        w1(i,j,k) = (1 - omega) * w_old(i,j,k) + omega * w_jacobi;
                    end
                end
            end

            % 检查收敛性
            if mod(iter, 50) == 0
                % 使用L2范数(均方根误差)来评估收敛性，它比最大值范数更稳定
                residual = sqrt(sum((w1(:) - w_old(:)).^2)) / sqrt(numel(w1));
                fprintf('迭代次数: %d, 残差 (L2 Norm): %e\n', iter, residual);

                if residual < tolerance
                    fprintf('收敛成功！\n');
                    break;
                end
            end
        end

        if iter == max_iter
            residual = sqrt(sum((w1(:) - w_old(:)).^2)) / sqrt(numel(w1));
            fprintf('达到最大迭代次数 %d，可能未完全收敛。最终残差: %e\n', max_iter, residual);
        end
        w(:,:,:,ii)=w1;
    end
        
end
