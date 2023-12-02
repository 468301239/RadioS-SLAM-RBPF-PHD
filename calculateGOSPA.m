function gospa = calculateGOSPA(trackerOutput, trueTargets, p, c)
    % trackerOutput: 跟踪器输出，每行包含一个目标的坐标 [x, y]
    % trueTargets: 真实目标的坐标，每行包含一个目标的坐标 [x, y]
    % p: GOSPA 参数，控制分布的稀疏性
    % c: GOSPA 参数，控制目标之间的相关性

    if isempty(trackerOutput) && isempty(trueTargets)
        gospa = 0;  % 两者都为空时，GOSPA为0
        return;
    end

    if isempty(trackerOutput) || isempty(trueTargets)
        gospa = inf;  % 其中一个为空时，GOSPA为无穷大
        return;
    end

    N = size(trackerOutput, 1);
    M = size(trueTargets, 1);

    % 计算距离矩阵
    distanceMatrix = pdist2(trackerOutput, trueTargets);

    % 计算匹配矩阵
    assignmentMatrix = zeros(N, M);
    for i = 1:N
        [~, minIndex] = min(distanceMatrix(i, :));
        assignmentMatrix(i, minIndex) = 1;
    end

    % 计算GOSPA
    miss = sum(1 - sum(assignmentMatrix, 2));  % 未匹配的目标数
    falseAlarm = sum(1 - sum(assignmentMatrix, 1));  % 多余的目标数
    mismatch = sum(sum(assignmentMatrix .* distanceMatrix));  % 匹配目标之间的距离

    gospa = p * (miss + falseAlarm) + c * mismatch;
end