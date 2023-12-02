function normalizedAngles = normalizeAngles(angles)
    % 输入：
    %   angles: 角度矩阵，可以是任意维度的矩阵
    % 输出：
    %   normalizedAngles: 标准化后的角度矩阵，范围在 -π 到 π 之间

    % 使用 mod 函数将角度限制在 0 到 2*pi 之间
    angles = mod(angles, 2*pi);

    % 将角度转换到 -pi 到 pi 之间
    normalizedAngles = angles;
    normalizedAngles(normalizedAngles > pi) = normalizedAngles(normalizedAngles > pi) - 2*pi;
end