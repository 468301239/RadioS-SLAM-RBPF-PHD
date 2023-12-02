clear;

%% hyper parameters
% landmark motion model
stateTransitionFcn = @constpos;
stateTransitionJacobianFcn=@constposjac;
hasAdditiveProcessNoise = true;
Q = 0.005;

% landmark measurement model
measurementFcn = @RngBrgMeasFcn;
measurementJacobianFcn = @RngBrgMeasFcnjac;
hasAdditiveMeasurementNoise = true;
Rngnoise=0.3;
Brgnoise=5*pi/180;

% birth strategy: sigma 2:0.9111, sigma 3: 0.9946
Rngthreshold=1;
Brgthreshold=20*pi/180;
wb=0.01;

% detection, survival and clutter intensity
Pd=0.9;
Ps=1;
Kc=0.000001;

% phd start
m0=zeros(2,0);
P0=zeros(2,2,0);
w0=zeros(1,0);

extractthreshold=0.09;
%% filter definition
phd = gmphd(m0,P0,...
    'Weights',w0,...
    'StateTransitionFcn',stateTransitionFcn,...
    'HasAdditiveProcessNoise',hasAdditiveProcessNoise,...
    'ProcessNoise',Q,...
    'MeasurementFcn',measurementFcn,...
    'HasAdditiveMeasurementNoise',hasAdditiveMeasurementNoise);

%% test
N=200;
history_phdstates=cell(N,1);
history_landmarks=cell(N,1);
history_phdStateCovariances=cell(N,1);
history_weights=cell(N,1);
load('data/meas.mat')
reservedMeasurementsnext=cell(0,1);
reservedlistsnext=zeros(1,0);
for iter=1:N
    measures=meas_complete{iter};
    predict(phd,1);
    
    % birth GCs    
    reservedMeasurements=reservedMeasurementsnext;
    reservedlists=reservedlistsnext;
    
    if ~isempty(reservedlists)
        birthPHD=getbirthPHD(reservedMeasurements,wb,diag([Rngnoise^2,Brgnoise^2]));
        append(phd, birthPHD);
    end
    
    if phd.NumComponents==0
        birthPHD=getbirthPHD(measures,wb,diag([Rngnoise^2,Brgnoise^2]));
        append(phd, birthPHD);
        history_phdstates{iter,1}=phd.States;
        history_landmarks{iter,1}=extractState(phd,extractthreshold);
        history_phdStateCovariances{iter,1}=phd.StateCovariances;
        history_weights{iter,1}=phd.Weights;
        continue;
    end
    
    reservedMeasurementsnext=cell(0,1);
    reservedlistsnext=zeros(1,0);
    allmeas=phd.MeasurementFcn(phd.States);
    for j=1:length(measures)
        measdeviation=(allmeas-measures{j}.Measurement);
        measdeviation(2,:)=normalizeAngles(measdeviation(2,:));
        measabsdeviation=abs(measdeviation);
        flagmatrix=[measabsdeviation(1,:)>Rngthreshold;measabsdeviation(2,:)>Brgthreshold];
        if  min(sum(flagmatrix,1))>0
            reservedMeasurementsnext=[reservedMeasurementsnext;measures(j)];
            reservedlistsnext=[j,reservedlistsnext];
        end
    end
    
    for j=1:length(reservedlistsnext)
        measures(reservedlistsnext(j))=[];
    end
    
    % start filtering process
    scale(phd, Ps);
    phdUndetected = clone(phd);
    scale(phdUndetected, 1 - Pd);
    if ~isempty(measures)
        phd.Detections = measures;
        detectionGroups = logical(eye(length(measures))); 
        logqij = likelihood(phd, detectionGroups);
        Lij = calculateScaling(logqij, phd.Weights, Pd, Kc);
        correct(phd, detectionGroups, Lij);
        append(phd, phdUndetected);
    else
        phd=phdUndetected;
    end    
    
    merge(phd,5);
    prune(phd, phd.Weights < 0.0001);
    
    history_phdstates{iter,1}=phd.States;
    history_landmarks{iter,1}=extractState(phd,extractthreshold);
    history_phdStateCovariances{iter,1}=phd.StateCovariances;
    history_weights{iter,1}=phd.Weights;
end


%% 画图
h=figure();
index=180;
states=history_phdstates{index};
statecovariances=history_phdStateCovariances{index};
weights=history_weights{index};
% 
% states=phd.States;
% statecovariances=phd.StateCovariances;
% weights=phd.Weights;

scatter(states(1,:),states(2,:))
hold on
scatter(target(1,:),target(2,:),'r+')
axis equal
box on;
grid on;
xlim([-50,50]);
ylim([-50,50]);
text(-20, 40, ['timeIndex=' num2str(index)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
[~,numcomponents]=size(states);
for i=1:numcomponents
%     3\sigma 置信度 0.9111
    plotGaussianEllipse(h, states(:,i), 4*statecovariances(:,:,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
end

for j=1:length(meas_complete{index})
    scatter(meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2)),meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2)),'g');
end

%% GOSPA 图，利用history_landmarks，即根据threshold被筛出来的GCs
tgm = trackGOSPAMetric('CutoffDistance',10,'Order',2,'Alpha',2,'SwitchingPenalty',0);
gospa=zeros(N,1);
for i=1:N
    track=zeros(2,0);
    for j=1:length(history_landmarks{i,1})
        track=[track,history_landmarks{i,1}(j).State];
    end
    [gospa(i),~,~]=GOSPA(target,track,2,5,2);
end
figure();
plot(gospa,'LineWidth',2,'Color',[1,0,0]);
box on;
grid on;
% g=figure();
% gridnum=30;
% % 生成x和y的均匀网格
% x = linspace(12.5, 13.5, gridnum);  % 100个点
% y = linspace(10.5, 11.5, gridnum);
% [x, y] = meshgrid(x, y);
% f=zeros(gridnum,gridnum);
% for i=1:gridnum
%     tic;
%     for j=1:gridnum
%         f(i,j)=calculatePHD(phd,[x(i,j);y(i,j)]);
%         title('三维曲面图');
%         xlabel('X轴');
%         ylabel('Y轴');
%         zlabel('f(x, y)');
%         colorbar;  % 显示颜色条
%     end
%     toc;
% end
% surf(x, y, f, 'EdgeColor', 'none');

%% 生成视频
% 步骤 2：设置视频参数
frameRate = 10;  % 帧率
duration = 20;    % 视频持续时间（秒）

% 步骤 3：创建视频对象
videoFile = 'output_video.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.FrameRate = frameRate;
open(videoObj);

% 步骤 4：循环生成每一帧
figure;
index=0;
h=figure();
for t = 0:1/frameRate:duration-1/frameRate
    % 使用 plot 函数绘制每一帧的内容
    clf(h);
    index=index+1;
    states=history_phdstates{index};
    statecovariances=history_phdStateCovariances{index};
    weights=history_weights{index};
    scatter(states(1,:),states(2,:))
    hold on
    scatter(target(1,:),target(2,:),'r+')
    axis equal
    box on;
    grid on;
    [~,numcomponents]=size(states);
    for i=1:numcomponents
    %     3\sigma 置信度 0.9111
        plotGaussianEllipse(h, states(:,i), 4*statecovariances(:,:,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
    end

    for j=1:length(meas_complete{index})
        scatter(meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2)),meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2)),'g');
    end
    scatter(0,0,'rs')
    title('Animated Plot');
    xlabel('X-axis');
    ylabel('Y-axis');

    text(-20, 40, ['timeIndex=' num2str(index)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
    % 将当前帧写入视频对象
    grid on;
    axis equal;
    xlim([-50,50]);
    ylim([-50,50]);
    
    writeVideo(videoObj, getframe(gcf));
end

% 步骤 5：保存视频并关闭对象
close(videoObj);
disp('Video generation complete.');

% 可选：播放生成的视频
% winopen('output_video.mp4');