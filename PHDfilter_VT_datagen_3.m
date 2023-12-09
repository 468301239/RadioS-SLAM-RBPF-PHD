clear all; close all; clc;
f=12.5;   % frequency hz
dt=1/f; % period s
T=30;   % exp time s
N=floor(T/dt);

fov=60;

bs=[0;0];
vt1=[0;10;0];
vt2=[10;-5;sqrt(125)];
vt3=[10;15;sqrt(125)];
vt4=[10;-5;sqrt(325)];

ax=0.5; ay=ax;    % acceleration for each direction m/s^2
pos_init=[0;0]; % initial position m
v_init=[1;0]/sqrt(1);   % initial velocity m/s

toarmse=0.3; % m
toalosrmse=0.05; % m
aoarmse=4*pi/180;   % rad
aoalosrmse=2*pi/180;
bias=0.3;
R=diag([toarmse^2,aoarmse^2]);

x=zeros(5,N);
x(:,1)=[pos_init;v_init;0.3];
A=eye(5);
A(1,3)=dt;
A(2,4)=dt;
B=[0.5*dt^2, 0;
    0, 0.5*dt^2;
    dt, 0;
    0, dt;
    0, 0];

% add dependencies
addpath('VT_datagen_dependencies');
%% groundtruth gen
for montecarlo=1:10
    flag=0;
    while (flag~=1)
        a=[normrnd(0,1,1,N)*ax;normrnd(0,1,1,N)*ay];
        for i=2:N
            x(:,i)=A*x(:,i-1)+B*a(:,i);
            if (x(2,i)>5)
                flag=0;
                break;
            end
            if(i==N)
                flag=1;
            end
        end
    end
    gtPose=[[0:0.08:29.92]',x'];
    %% plot
    h=figure();
    grid on;
    hold on;
    xlim([-5,30]);
    ylim([-10,20]);
    scatter([bs(1);vt1(1);vt2(1);vt3(1);vt4(1)],[bs(2);vt1(2);vt2(2);vt3(2);vt4(2)]);
    plot(x(1,:),x(2,:));

    %% measurement gen
    lambda = 0.02; % expectation of the Poisson distribution
    P_D=0.95;      % probability of detection

    bsmeas=zeros(2,N);  % 先toa后dod最后aoa
    vt1meas=zeros(3,N);
    vt2meas=zeros(3,N);
    vt3meas=zeros(3,N);
    vt4meas=zeros(3,N);

    
    bsmeascpp=zeros(2,N);  % 先toa后aoa
    vt1meascpp=zeros(2,N);
    vt2meascpp=zeros(2,N);
    vt3meascpp=zeros(2,N);
    vt4meascpp=zeros(2,N);
    
    bsmeasoriginal=zeros(3,N);  % 先toa后dod最后aoa
    vt1measoriginal=zeros(3,N);
    vt2measoriginal=zeros(3,N);
    vt3measoriginal=zeros(3,N);
    vt4measoriginal=zeros(3,N);

    bsmeascpporiginal=zeros(2,N);  % 先toa后aoa
    vt1meascpporiginal=zeros(2,N);
    vt2meascpporiginal=zeros(2,N);
    vt3meascpporiginal=zeros(2,N);
    vt4meascpporiginal=zeros(2,N);

    gtPosetrans=gtPose';
    for i=1:N
        alpha=GetAngle(x(3:4,i),[0;0]);
        bsmeasoriginal(:,i)=[norm(x(1:2,i)-bs)+bias;getAOD([0;0;0],x(:,i),[0;0;0],'LOS');getAOA([0;0;0],x(:,i),[0;0;0],'LOS',alpha)'];
        vt1measoriginal(:,i)=[norm(x(1:2,i)-vt1(1:2))+vt1(3)+bias;getAOD([0;0;0],x(:,i),vt1,'VA');getAOA([0;0;0],x(:,i),vt1,'VA',alpha)'];
        vt2measoriginal(:,i)=[norm(x(1:2,i)-vt2(1:2))+vt2(3)+bias;getAOD([0;0;0],x(:,i),vt2,'VA')';getAOA([0;0;0],x(:,i),vt2,'VA',alpha)'];
        vt3measoriginal(:,i)=[norm(x(1:2,i)-vt3(1:2))+vt3(3)+bias;getAOD([0;0;0],x(:,i),vt2,'SP')';getAOA([0;0;0],x(:,i),vt3,'SP',alpha)'];
        vt4measoriginal(:,i)=[norm(x(1:2,i)-vt4(1:2))+vt4(3)+bias;getAOD([0;0;0],x(:,i),vt3,'VA')';getAOA([0;0;0],x(:,i),vt4,'VA',alpha)'];
        
        bsmeascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-bs)+bias;GetAngle([0,0],gtPosetrans(2:3,i))];
        vt1meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt1(1:2))+vt1(3)+bias;GetAngle(vt1(1:2),gtPosetrans(2:3,i))];
        vt2meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt2(1:2))+vt2(3)+bias;GetAngle(vt2(1:2),gtPosetrans(2:3,i))];
        vt3meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt3(1:2))+vt3(3)+bias;GetAngle(vt3(1:2),gtPosetrans(2:3,i))];
        vt4meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt4(1:2))+vt4(3)+bias;GetAngle(vt4(1:2),gtPosetrans(2:3,i))];
    end
    
    for count=1:10
        % white noise addup
        bsnoise=[randn(1,N)*toalosrmse;randn(1,N)*aoalosrmse;randn(1,N)*aoalosrmse];
        vt1noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt2noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt3noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt4noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];

        bsmeas=bsmeasoriginal+bsnoise;
        vt1meas=vt1measoriginal+vt1noise;
        vt2meas=vt2measoriginal+vt2noise;
        vt3meas=vt3measoriginal+vt3noise;
        vt4meas=vt4measoriginal+vt4noise;
        
        bsmeascpp=bsmeascpporiginal+bsnoise([1,3],:);
        vt1meascpp=vt1meascpporiginal+vt1noise([1,3],:);
        vt2meascpp=vt2meascpporiginal+vt2noise([1,3],:);
        vt3meascpp=vt3meascpporiginal+vt3noise([1,3],:);
        vt4meascpp=vt4meascpporiginal+vt4noise([1,3],:);

        bsmeascpp(2,:)=normalizeAngles(bsmeascpp(2,:));
        vt1meascpp(2,:)=normalizeAngles(vt1meascpp(2,:));
        vt2meascpp(2,:)=normalizeAngles(vt2meascpp(2,:));
        vt3meascpp(2,:)=normalizeAngles(vt3meascpp(2,:));
        vt4meascpp(2,:)=normalizeAngles(vt4meascpp(2,:));

        bsmeascpp(1,:)=normalizeRanges(bsmeascpp(1,:));
        vt1meascpp(1,:)=normalizeRanges(vt1meascpp(1,:));
        vt2meascpp(1,:)=normalizeRanges(vt2meascpp(1,:));
        vt3meascpp(1,:)=normalizeRanges(vt3meascpp(1,:));
        vt4meascpp(1,:)=normalizeRanges(vt4meascpp(1,:));

        h1=figure();
        hold on;
        plot(vt1meascpp(1,:));
        plot(vt2meascpp(1,:));
        plot(vt3meascpp(1,:));
        plot(vt4meascpp(1,:));
        %% code from GenData
        % hyper parameter.
        parameter.N=N;
        parameter.T=T;
        parameter.dT=parameter.T/parameter.N;

        sigma.DOA_az = aoarmse; 
        sigma.DOA_el = 0.00001; 
        sigma.DOD_az = aoarmse; 
        sigma.DOD_el = 0.00001; 
        sigma.TOA = toarmse; 
        % [TOA, DODaz, DODel, DOAaz, DOAel]
        measurementCovariance = diag([sigma.TOA^2,sigma.DOD_az^2,sigma.DOD_el^2,sigma.DOA_az^2,sigma.DOA_el^2]);

        % 为了适应matlab代码格式
        % para
        para.save = 1; % 1) save;
        para.MC = 10; % # Monte Carlo run
        para.TIME = parameter.N; % # time evolution (para.TIME = 2*pi/Vel_rot.ini/para.t_du)
        para.t_du = parameter.dT; % time duration
        para.N_vehicle = 1; % number of vehicles
        para.N_VA = 1; % # number of VAs
        para.N_SP = 3;  % # number of SPs
        para.SPVisibilityRadius=fov; % FoV for SP, [m]
        para.Rmax = 100; % maximum distance range, [m]
        
        para.lambda = lambda; % expectation of the Poisson distribution
        para.P_D=P_D;      % probability of detection
        
        
        % units: [x m; y m; z m; theta rad; v m/s; theta dot rad/s; B m]^2 

        % process noise: to be modified

        para.ax=ax;
        para.ay=ay;
        para.ProcessNoiseCovariance=[];     % inaccurate expression
        para.InitialState = [x(1:2,1)',0,GetAngle(x(3:4,1),[0;0]),norm(x(3:4,1)),0,0];

        % Channel
        Channel.Values = zeros(5,para.N_VA+para.N_SP+1,para.TIME); % ordered: TOA, DOD(az,el) DOA(az,el)
        Channel.Labels=['TOA ' ' DODaz' ' DODel' ' DOAaz' ' DOAel'];
        Channel.Visible.TOT= ones(para.N_VA+para.N_SP+1,para.TIME,para.N_vehicle);  % 0 or 1
        Channel.Visible.LOS= ones(1,para.TIME,para.N_vehicle);  % 0 or 1
        Channel.Visible.VA= ones(para.N_VA,para.TIME,para.N_vehicle);  % 0 or 1
        Channel.Visible.SP= ones(para.N_SP,para.TIME,para.N_vehicle);  % 0 or 1
        Channel.Clutter=[];

        % state translate
        % units: [x m; y m; z m; theta rad; v m/s; theta dot rad/s; B m]
        state=[x(1:2,:);zeros(1,parameter.N);zeros(1,parameter.N);sqrt(x(3,:).^2+x(4,:).^2);zeros(1,parameter.N);zeros(1,parameter.N)];
        for i=1:para.TIME
            state(4,i)=GetAngle(x(3:4,i),[0;0]);
        end

        % [TOA, DODaz, DODel, DOAaz, DOAel]
        for i=1:para.TIME
            Channel.Values(:,1,i)=[bsmeas(1,i);bsmeas(2,i);0;bsmeas(3,i);0];

            Channel.Values(:,2,i)=[vt1meas(1,i);vt1meas(2,i);0;vt1meas(3,i);0];

            Channel.Values(:,3,i)=[vt2meas(1,i);vt2meas(2,i);0;vt2meas(3,i);0];
            Channel.Values(:,4,i)=[vt3meas(1,i);vt3meas(2,i);0;vt3meas(3,i);0];
            Channel.Values(:,5,i)=[vt4meas(1,i);vt4meas(2,i);0;vt4meas(3,i);0];
            countall=-1;
            if i<=para.TIME/5
                v.Time(i).measurement(:,1)=Channel.Values(:,1,i);
                v.Time(i).measurement(:,2)=Channel.Values(:,2,i);
                v.Time(i).measurement(:,3)=Channel.Values(:,3,i);
                v.Time(i).measurement(:,4)=Channel.Values(:,4,i);
                v.Time(i).measurement(:,5)=Channel.Values(:,5,i);
                countall=5;
            else
                Channel.Visible.LOS(1,i)=0;
                Channel.Visible.TOT(1,i)=0;
                v.Time(i).measurement(:,1)=Channel.Values(:,2,i);
                v.Time(i).measurement(:,2)=Channel.Values(:,3,i);
                v.Time(i).measurement(:,3)=Channel.Values(:,4,i);
                v.Time(i).measurement(:,4)=Channel.Values(:,5,i);
                countall=4;
            end

            N_clutter = poissrnd(para.lambda,1);
            Channel.Clutter(i).vehicle(1).visible = N_clutter;
            
            if N_clutter > 0
                for m = 1:N_clutter
                    countall=countall+1;
                    Channel.Clutter(i).vehicle(1).measurement(:,m) = [para.Rmax*rand; 2*pi*rand; 0; 2*pi*rand; 0];
                    cluttertemp(i).measurement(:,m)=[Channel.Clutter(i).vehicle(1).measurement(1,m),angcali(Channel.Clutter(i).vehicle(1).measurement(4,m)+state(4,i))];
                    v.Time(i).measurement(:,countall)=Channel.Clutter(i).vehicle(1).measurement(:,m);
                end
            end

            v.Time(i).measurement=transpose(v.Time(i).measurement);
        end
        % units: [x m; y m; z m; theta rad; v m/s; theta dot rad/s; B m]^2 
        BS.pos=[0;0;0];
        VA(1).pos=[0; 10; 0];
        SP(1).pos=[10; -5; 0];
        SP(2).pos=[10; 15; 0];

        if para.save ==1
            mkdir(['FinalData_4VT2/Meas_',num2str(montecarlo)]);
            Sp = ['FinalData_4VT2/Meas_',num2str(montecarlo),'/measurement_',num2str(count),'_' num2str(para.TIME)];
            filepath = Sp;
            save(filepath, 'state', 'para', 'BS', 'VA', 'SP', 'Channel', 'v', 'measurementCovariance')
            
            vtmeasurements=cell({vt1meascpp,vt2meascpp,vt3meascpp,vt4meascpp});
            % saving
            meas_complete=cell(N,1);

            N_clutter = poissrnd(lambda,1);

            for i=1:N
                meas_complete{i}=cell(0,1);
                mp = struct(OriginPosition = [gtPose(i,2:3), bias]);
                for j=1:4
                    if rand()<P_D && vtmeasurements{j}(1,i)<fov                
                        tempcell=cell(1,1);
                        tempcell{1,1}=objectDetection((i-1)*dt,vtmeasurements{j}(:,i),'MeasurementNoise',R,MeasurementParameters=mp);
                        meas_complete{i}=[meas_complete{i};tempcell];
                    end
                end

                for j=1:N_clutter
                    r_clutter=fov*rand();
                    b_clutter=2*pi*rand()-pi;
                    tempcell=cell(1,1);
                    tempcell{1,1}=objectDetection((i-1)*dt,[r_clutter;b_clutter],'MeasurementNoise',R,MeasurementParameters=mp);
                    meas_complete{i}=[meas_complete{i};tempcell];
                end
            end

            target=[vt1,vt2,vt3,vt4];
            save(['FinalData_4VT2/Meas_',num2str(montecarlo),'/meas_',num2str(count),'_phd.mat'],'meas_complete','target','gtPose','bsmeascpp');
            print (h1, '-dpng', ['FinalData_4VT2/Meas_',num2str(montecarlo),'/meas_',num2str(count),'.png']);
            close(h1);
        end
        print (h, '-dpng', ['FinalData_4VT2/Meas_',num2str(montecarlo),'/gtTraj.png']);
        save( ['FinalData_4VT2/Meas_' num2str(montecarlo) '/gtPose.mat'],'gtPose');
        clear parameter;
        clear Channel;
        clear v;
        
    end
end

