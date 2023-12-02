clear;
%% filtering
filter=myphd();

N=375;
history_phdstates=cell(N,1);
history_filters=cell(N,1);
history_phdStateCovariances=cell(N,1);
history_weights=cell(N,1);
history_landmarks=cell(N,1);



load('dataVT/meas.mat')
% load('test.mat');
% for iter=1:N
%     for count=1:length(meas_complete{iter})
%         meas_complete{iter}{count}.MeasurementParameters.OriginPosition=[estimatedStates(1,iter),estimatedStates(2,iter),estimatedStates(5,iter)];
%     end
% end

for iter=1:N
    measures=meas_complete{iter};
    filter=filter.prdupd(measures);
    
    history_filters{iter}=clone(filter);
    
    history_phdstates{iter,1}=filter.phd.States;
    history_phdStateCovariances{iter,1}=filter.phd.StateCovariances;
    history_weights{iter,1}=filter.phd.Weights;
    history_landmarks{iter,1}=extractState(filter.phd,filter.config.extractthreshold);

end

ending=0;

%% 画图
h=figure();
index=38;
states=history_filters{index}.phd.States;
statecovariances=history_filters{index}.phd.StateCovariances;
weights=history_filters{index}.phd.Weights;

% states=phd.States;
% statecovariances=phd.StateCovariances;
% weights=phd.Weights;

scatter(states(1,:),states(2,:))
hold on
scatter(target(1,:),target(2,:),'r+')
axis equal
box on;
grid on;
xlim([-20,25]);
ylim([-25,20]);
text(-20, 40, ['Step=' num2str(index)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
[~,numcomponents]=size(states);
for i=1:numcomponents
%     3\sigma 置信度 0.9111
    plotGaussianEllipse(h, states(1:2,i), 4*statecovariances(1:2,1:2,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
end

for j=1:length(meas_complete{index})
    scatter(meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2)),meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2)),'g');
end

%% GOSPA 图，利用history_landmarks，即根据threshold被筛出来的GCs
gospa=zeros(N,1);
for i=1:N
    track=zeros(2,0);
    for j=1:length(history_landmarks{i,1})
        track=[track,history_landmarks{i,1}(j).State(1:2)];
    end
    [gospa(i),~,~]=GOSPA(target(1:2,:),track,2,5,2);
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
duration = N/frameRate;    % 视频持续时间（秒）

% 步骤 3：创建视频对象
videoFile = 'output_video.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
videoObj.FrameRate = frameRate;
open(videoObj);

% 步骤 4：循环生成每一帧
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

    scatter([0],[0],70,'filled','b');
%     scatter([10],[-5],70,'filled','rs');
    hold on
    box on;
    grid on;
    [~,numcomponents]=size(states);
    for i=1:numcomponents
    %     2\sigma 置信度 0.8696
        plotGaussianEllipse(h, states(1:2,i), 4*statecovariances(1:2,1:2,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
    end

    for j=1:length(meas_complete{index})
        h1=plot([gtPose(index,2),gtPose(index,2)+meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2))],[gtPose(index,3),gtPose(index,3)+meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2))],'Color',[0,1,0],'LineWidth',2);
    end
%     scatter(0,0,'rs')
    
    
    hold on;
    plot([-1000,10,10],[5,5,-1000],'color',[0,0,0],'linewidth',2);

    theta=30*pi/180;
    d=2;
    for x=-25:1.5:10
        plot([x,x+d*cos(theta)],[5,5+d*sin(theta)],'color',[0,0,0],'linewidth',1);
    end
    for y=-25:1.5:5
        plot([10,10+d*cos(theta)],[y,y+d*sin(theta)],'color',[0,0,0],'linewidth',1);
    end

    h3=plot(gtPose(1:index,2),gtPose(1:index,3),'color',[0,0,1],'linewidth',1.5);
    h4=plot(estimatedStates(1,1:index),estimatedStates(2,1:index),'color',[1,0,0],'linewidth',1.5);
    scatter(estimatedStates(1,index),estimatedStates(2,index),'filled','r');
    
%     vts=scatter([0,10],[10,15],70,[1, 0, 0],'filled','^');
%     vts.MarkerFaceAlpha=0.5;
%     vts.MarkerEdgeAlpha=0.5;

%     text(41,3,"Wall",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','right','Interpreter','latex');
%     text(0-1.5,-2,"Anchor",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
%     text(0-1.5,10-2,"VT1",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
%     text(10-7.5,-5-2,"Scatterer (VT2, VT4)",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
%     text(10-1.5,15-2,"VT3",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
    text(-10, 20, ['Step:' num2str(index) ' Pd: 0.95'], 'FontSize', 10, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'Color', [1, 0, 0]);
    text(-10, 17, ['Ellipse: 2 sigma'], 'FontSize', 10, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'right', 'Interpreter', 'latex', 'Color', [1, 0, 0]);
    set(gca,'FontName','Times New Roman');

    scatter(target(1,:),target(2,:),'r+')
    grid on;
    box on;
    xlabel("x (m)",'Interpreter','latex','FontSize',16);
    ylabel("y (m)",'Interpreter','latex','FontSize',16);
    
    title("PHD SLAM")
    legend([h1,h3,h4],"Measurements","Trajectory (GT)","Trajectory (Estimated)",'Location', 'southwest');

    xlim([-25,30]);
    ylim([-25,30]);
    set(gcf,'unit','centimeters','position',[1,2,16,16*3/4]);
    set(gca,'FontSize',16)
    set(gca,'TickLabelInterpreter','latex');
   
    
    writeVideo(videoObj, getframe(gcf));
end

% 步骤 5：保存视频并关闭对象
close(videoObj);
disp('Video generation complete.');

% 可选：播放生成的视频
% winopen('output_video.mp4');
%%

figure();
hold on;
plot([-1000,1000],[5,5],'color',[0,0,0],'linewidth',2);

theta=30*pi/180;
d=2;
for x=-10:1.5:45
    plot([x,x+d*cos(theta)],[5,5+d*sin(theta)],'color',[0,0,0],'linewidth',1.5);
end

h3=plot(gtPose(:,2),gtPose(:,3),'color',[0,0,1],'linewidth',2);

% vts=scatter([0,10],[10,15],70,[1, 0, 0],'filled','^');
% vts.MarkerFaceAlpha=0.5;
% vts.MarkerEdgeAlpha=0.5;

scatter([0],[0],70,'filled','b');
scatter([10],[-5],70,'filled','rs');
text(41,3,"Wall",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','right','Interpreter','latex');
text(0-1.5,-2,"Anchor",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
text(0-1.5,10-2,"VT1",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
text(10-7.5,-5-2,"Scatterer (VT2, VT4)",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
text(10-1.5,15-2,"VT3",'FontSize',14,'FontName','Times New Roman','HorizontalAlignment','left','Interpreter','latex');
set(gca,'FontName','Times New Roman');

grid on;
box on;
xlabel("x (m)",'Interpreter','latex','FontSize',16);
ylabel("y (m)",'Interpreter','latex','FontSize',16);
title("PHD Mapping with GT Trajectory")
legend([h3],"Trajectory (GT)");

xlim([-3,43]);
ylim([-13,18]);
set(gcf,'unit','centimeters','position',[1,2,16,16*3/4]);
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex');