function [dsdz] = model_gradient_z(s,z)
    % 单独计算垂向梯度
    %
    % 输入：
    % s - 变量
    % z - z坐标,一维
    %
    % 输出：
    % dsdz - z方向梯度

    % 地球半径（单位：米）
    dsdz=zeros(size(s));
    % 中间点 使用中心差分
    for j=2:size(s,3)-1
        dsdz(:,:,j)=(s(:,:,j+1)-s(:,:,j-1)) ./ (z(j+1)-z(j-1));
    end

    % 边界点 使用单侧差分
    dsdz(:,:,1)=(s(:,:,2)-s(:,:,1)) ./ (z(2)-z(1));
    dsdz(:,:,end)=(s(:,:,end)-s(:,:,end-1)) ./ (z(end)-z(end-1));

end