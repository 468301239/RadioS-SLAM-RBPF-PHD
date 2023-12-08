clear all;
close all;
clc;

% myfilter=history_filters{375};
%% parameters
N_p=100;
state_init=[8;-10;0;1;0.3];
weight_init=1/N_p;
N_eff_threshold=50;

% load measurements
for trajcount=1:2
    for meascount=1:10
        load(['FinalData/Meas_',num2str(trajcount),'/meas_',num2str(meascount),'_phd.mat'])
        % load('test.mat');
        % priorstates=estimatedStates;
        N=length(meas_complete);
        estimatedStates=zeros(5,N);
        history_phdstates=cell(N,1);
        history_filters=cell(N,1);
        history_phdStateCovariances=cell(N,1);
        history_weights=cell(N,1);
        history_landmarks=cell(N,1);

        % initialize particles
        iter=1;
        ParticleSet = cell(N_p,1);
        for i = 1:N_p
            ParticleSet{i} = myparticle(state_init, weight_init);
        end
        measures=meas_complete{1};
        for i=1:N_p
            ParticleSet{i}=setmeas(ParticleSet{i},measures);
            ParticleSet{i}=mapiter(ParticleSet{i});
        %     ParticleSet{i}.map=clone(myfilter);
        end
        estimatedStates(:,1)=[8;-10;0;1;0.3];
        history_filters{iter}=ParticleSet{1}.map;
        history_phdstates{iter,1}=ParticleSet{1}.map.phd.States;
        history_phdStateCovariances{iter,1}=ParticleSet{1}.map.phd.StateCovariances;
        history_weights{iter,1}=ParticleSet{1}.map.phd.Weights;
        history_landmarks{iter,1}=extractState(ParticleSet{1}.map.phd,ParticleSet{1}.map.config.extractthreshold);

        for iter=2:N
            % particle propagation
            for i=1:N_p
                ParticleSet{i}=predict(ParticleSet{i});
        %         if iter<=N
        %             ParticleSet{i}.state=priorstates(:,iter);
        %         end
            end

            % map predict and update
            measures=meas_complete{iter};
            for i=1:N_p
                ParticleSet{i}=setmeas(ParticleSet{i},measures);
                ParticleSet{i}=mapiter(ParticleSet{i});
        %         ParticleSet{i}.map=clone(myfilter);
            end

            % importance weighting
            weightsum=0;
            for i=1:N_p
                ParticleSet{i}=reweight(ParticleSet{i});
        %         xtemp=ParticleSet{i}.state(1);
        %         ytemp=ParticleSet{i}.state(2);
        %         ParticleSet{i}.weight=ParticleSet{i}.weight*multifeaturereweight(ParticleSet{i}.measures,ParticleSet{i}.map,[xtemp,ytemp,0.3],100);
                if iter<=N/5
                    ParticleSet{i}=reweightlos(ParticleSet{i},bsmeascpp(:,iter));
        %             z=bsmeascpp(:,iter);
        %             mp = struct(OriginPosition = ParticleSet{i}.state');
        %             zexp=RngBrgMeasFcnVT([0;0;0],mp);
        %             zCov=diag([ParticleSet{i}.map.config.Rngnoiselos^2,ParticleSet{i}.map.config.Brgnoiselos^2]);
        %             error = z-zexp;
        %             error(2)=normalizeAngles(error(2));
        %             md2=error'/zCov*error;
        %             ParticleSet{i}.weight=ParticleSet{i}.weight*1/sqrt(det(2*pi*zCov))*exp(-1/2*md2);
                end
                weightsum=weightsum+ParticleSet{i}.weight;
            end
            for i=1:N_p
                ParticleSet{i}.weight=ParticleSet{i}.weight/weightsum;
            end

            % calculate effective particle number and resample
            weightlist=[];
            for i=1:N_p
                weightlist=[weightlist,ParticleSet{i}.weight];
            end

            % extract estimated state and map
            maxdis=0;
            maxindex=1;
            for i=1:N_p
        %         dis=norm(ParticleSet{i}.state-gtPose(iter,2:6)');
        %         if dis>maxdis
        %             maxindex=i;
        %             maxdis=dis;
        %         end
                estimatedStates(:,iter)=estimatedStates(:,iter)+ParticleSet{i}.weight*ParticleSet{i}.state;
            end
        %     estimatedStates(:,iter)=ParticleSet{maxindex}.state;

            [~,maxindex]=max(weightlist);
            history_filters{iter}=clone(ParticleSet{maxindex}.map);
            history_phdstates{iter,1}=ParticleSet{maxindex}.map.phd.States;
            history_phdStateCovariances{iter,1}=ParticleSet{maxindex}.map.phd.StateCovariances;
            history_weights{iter,1}=ParticleSet{maxindex}.map.phd.Weights;
            history_landmarks{iter,1}=extractState(ParticleSet{maxindex}.map.phd,ParticleSet{maxindex}.map.config.extractthreshold);

            N_eff=1/(weightlist*weightlist');
            if N_eff<N_eff_threshold
                cumulative_weights = cumsum(weightlist);
                random_numbers = rand(1, N_p);
                ParticleSet_resampled = cell(N_p, 1);
                for i = 1:N_p
                    selected_index = find(cumulative_weights >= random_numbers(i), 1);
                    ParticleSet_resampled{i} = clone(ParticleSet{selected_index});
                    ParticleSet_resampled{i}.weight=1/N_p;
                end
                ParticleSet=cell(N_p,1);
                for i = 1:N_p
                    ParticleSet{i}=clone(ParticleSet_resampled{i});
                end
            end
            fprintf('Time: %s/375\n', num2str(iter)); 
        end
        
        % 画轨迹图
        h=figure();
        plot(estimatedStates(1,:),estimatedStates(2,:));
        hold on;
        plot(gtPose(:,2),gtPose(:,3));
        grid on;
        box on;
        axis equal;
        estPose=estimatedStates';
        serror=zeros(N,1);
        for i=1:N
            serror(i)=(estPose(i,1:2)-gtPose(i,2:3))*(estPose(i,1:2)-gtPose(i,2:3))';
        end
        rmse=sqrt(mean(serror));

        % 画map图
        h1=figure();
        index=N;
        states=history_phdstates{index};
        statecovariances=history_phdStateCovariances{index};
        weights=history_weights{index};

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
            plotGaussianEllipse(h1, states(1:2,i), 4*statecovariances(1:2,1:2,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
        end

        for j=1:length(meas_complete{index})
            scatter(meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2)),meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2)),'g');
        end

        mkdir(['FinalResult/Meas_',num2str(trajcount)]);
        print (h, '-dpng', ['FinalResult/Meas_',num2str(trajcount),'/result_',num2str(meascount),'_phd_traj.jpg']);
        print (h1, '-dpng', ['FinalResult/Meas_',num2str(trajcount),'/result_',num2str(meascount),'_phd_map.jpg']);
        close(h);
        close(h1);

        Sp = ['FinalResult/Meas_',num2str(trajcount),'/result_',num2str(meascount),'_phd.mat'];
        filepath = Sp;
        save(filepath, 'rmse','estimatedStates', 'history_phdstates', 'history_filters', 'history_phdStateCovariances', 'history_weights', 'history_landmarks');
    end
end

%     %% 画map图
%     h=figure();
%     index=iter-1;
%     states=history_phdstates{index};
%     statecovariances=history_phdStateCovariances{index};
%     weights=history_weights{index};
% 
%     % states=phd.States;
%     % statecovariances=phd.StateCovariances;
%     % weights=phd.Weights;
% 
%     scatter(states(1,:),states(2,:))
%     hold on
%     scatter(target(1,:),target(2,:),'r+')
%     axis equal
%     box on;
%     grid on;
%     xlim([-20,25]);
%     ylim([-25,20]);
%     text(-20, 40, ['Step=' num2str(index)], 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'red');
%     [~,numcomponents]=size(states);
%     for i=1:numcomponents
%     %     3\sigma 置信度 0.9111
%         plotGaussianEllipse(h, states(1:2,i), 4*statecovariances(1:2,1:2,i), 'MeanColor', 'b', 'EllipseColor', 'm', 'LineWidth', 1, 'Fill', true, 'Alpha', min(0.3*weights(i),1));
%     end
% 
%     for j=1:length(meas_complete{index})
%         scatter(meas_complete{index}{j}.Measurement(1)*cos(meas_complete{index}{j}.Measurement(2)),meas_complete{index}{j}.Measurement(1)*sin(meas_complete{index}{j}.Measurement(2)),'g');
%     end
%     %% 画轨迹图
%     plot(estimatedStates(1,:),estimatedStates(2,:));
%     hold on;
%     plot(gtPose(:,2),gtPose(:,3));
%     grid on;
%     box on;
%     axis equal;
%     estPose=estimatedStates';
%     serror=zeros(375,1);
%     for i=1:375
%         serror(i)=(estPose(i,1:2)-gtPose(i,2:3))*(estPose(i,1:2)-gtPose(i,2:3))';
%     end
%     rmse=sqrt(mean(serror));
% 
% %% test
% index=10;
% states=history_phdstates{index};
% weightlist=history_weights{index};
% % for i=1:N_p
% %     states=[states,ParticleSet{i}.state];
% %     weightlist=[weightlist,ParticleSet{i}.weight];
% % end
% 
% colors = weightlist / max(weightlist);
% 
% % 可视化
% scatter(states(1, :), states(2, :), 50, weightlist, 'filled');
% colorbar; % 添加颜色条
% title('粒子分布可视化');
% xlabel('X轴');
% ylabel('Y轴');