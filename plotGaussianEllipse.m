function plotGaussianEllipse(h, x, P, varargin)
    % 1D: （μ-σ,μ+σ)中的概率为0.6827；数值分布在（μ-2σ,μ+2σ)中的概率为0.9545；数值分布在（μ-3σ,μ+3σ)中的概率为0.9973
    % 2D: sigma 1: 0.466, sigma 2:0.9111, sigma 3: 0.9946
    % 3D: sigma 1: 0.318, sigma 2:0.8696, sigma 3: 0.9919
    % 生成椭圆数据
    numPoints = 100;
    theta = linspace(0, 2*pi, numPoints);
    points = [cos(theta); sin(theta)];

    % 计算椭圆的半长轴和半短轴
    [U, S, ~] = svd(P);
    lambda = sqrt(diag(S));

    % 将椭圆数据旋转并缩放
    ellipsePoints = U * (lambda .* points);

    % 平移椭圆到指定均值
    ellipsePoints = bsxfun(@plus, ellipsePoints, x);

    % 解析可选参数
    p = inputParser;
    addParameter(p, 'MeanColor', 'ro');  % 默认均值点颜色是红色
    addParameter(p, 'EllipseColor', 'b');  % 默认椭圆颜色是蓝色
    addParameter(p, 'LineWidth', 1.5);  % 默认线宽是1.5
    addParameter(p, 'Fill', false);  % 默认不填充椭圆
    addParameter(p, 'Alpha', 1);  % 默认不透明

    parse(p, varargin{:});

    % 绘制椭圆
    figure(h);
    hold on;

    if p.Results.Fill
        fill(ellipsePoints(1, :), ellipsePoints(2, :), p.Results.EllipseColor, 'LineWidth', p.Results.LineWidth, 'FaceAlpha', p.Results.Alpha);
    else
        plot(x(1), x(2), p.Results.MeanColor);  % 绘制均值点
        plot(ellipsePoints(1, :), ellipsePoints(2, :), 'Color', p.Results.EllipseColor, 'LineWidth', p.Results.LineWidth);
    end

%     title('高斯分布随机变量的椭圆可视化');
%     xlabel('X轴');
%     ylabel('Y轴');
end