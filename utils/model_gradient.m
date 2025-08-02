function [dsdx,dsdy] = model_gradient(x,y,s)
    % 计算不均匀网格下的梯度
    % 因为模式网格文件的坐标不均匀而开发
    %
    % 输入：
    % x - 网格文件中读取的x坐标，注意单位m (二维数组) 算不出来就是没转置
    % y - 网格文件中读取的x坐标，注意单位m (二维数组) 算不出来就是没转置
    %
    % 输出：
    % dsdx - x方向梯度
    % dsdy - y方向梯度

    % 地球半径（单位：米）
    R = 6371000;
    dsdx=zeros(size(s)); dsdy=zeros(size(s));
    % 中间点 使用中心差分
    hnum=size(s,3);
    for j=2:size(x,2)-1
        dsdy(:,j,:)=(s(:,j+1,:)-s(:,j-1,:)) ./ repmat((y(:,j+1)-y(:,j-1)),[1 1 hnum]);
    end
    for i=2:size(x,1)-1
        dsdx(i,:,:)=(s(i+1,:,:)-s(i-1,:,:)) ./ repmat((x(i+1,:)-x(i-1,:)),[1 1 hnum]);
    end
    
    % 边界点 使用单侧差分
    dsdy(:,1,:)=(s(:,2,:)-s(:,1,:)) ./ repmat((y(:,2)-y(:,1)),[1 1 hnum]);
    dsdy(:,end,:)=(s(:,end,:)-s(:,end-1,:)) ./ repmat((y(:,end)-y(:,end-1)),[1 1 hnum]);
    
    dsdx(1,:,:)=(s(2,:,:)-s(1,:,:)) ./ repmat((x(2,:)-x(1,:)),[1 1 hnum]);
    dsdx(end,:,:)=(s(end,:,:)-s(end-1,:,:)) ./ repmat((x(end,:)-x(end-1,:)),[1 1 hnum]);
end