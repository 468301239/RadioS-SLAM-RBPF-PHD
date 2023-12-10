clear all; close all; clc;
f=12.5;   % frequency hz
dt=1/f; % period s
T=30;   % exp time s
N=floor(T/dt);

fov=30;

bs=[0;0];
vt1=[0;10;0];
vt2=[20;0;0];
vt3=[5;-5;5*sqrt(2)];
vt4=[20;10;0];
vt5=[5;15;5*sqrt(2)];
vt6=[15;-5;5*sqrt(2)];
vt7=[5;-5;5*sqrt(10)];
vt8=[5;-5;5*sqrt(10)];

ax=0.5; ay=ax;    % acceleration for each direction m/s^2
pos_init=[8;-10]; % initial position m
v_init=[0;1]/sqrt(1);   % initial velocity m/s

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
for montecarlo=3:10
    flag=0;
    while (flag~=1)
        a=[normrnd(0,1,1,N)*ax;normrnd(0,1,1,N)*ay];
        for i=2:N
            x(:,i)=A*x(:,i-1)+B*a(:,i);
            if or(x(2,i)>5,x(1,i)>10)
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
    xlim([-20,25]);
    ylim([-25,20]);
    scatter([bs(1);vt1(1);vt2(1);vt3(1);vt4(1);vt5(1);vt6(1);vt7(1);vt8(1)],...
    [bs(2);vt1(2);vt2(2);vt3(2);vt4(2);vt5(2);vt6(2);vt7(2);vt8(2)]);
    plot(x(1,:),x(2,:),'r');
    plot([-100,10,10],[5,5,-100],'b');

    %% measurement gen
    lambda = 0.02; % expectation of the Poisson distribution
    P_D=0.95;      % probability of detection

    bsmeas=zeros(2,N);  % 先toa后dod最后aoa
    vt1meas=zeros(3,N);
    vt2meas=zeros(3,N);
    vt3meas=zeros(3,N);
    vt4meas=zeros(3,N);
    vt5meas=zeros(3,N);
    vt6meas=zeros(3,N);
    vt7meas=zeros(3,N);
    vt8meas=zeros(3,N);

    
    bsmeascpp=zeros(2,N);  % 先toa后aoa
    vt1meascpp=zeros(2,N);
    vt2meascpp=zeros(2,N);
    vt3meascpp=zeros(2,N);
    vt4meascpp=zeros(2,N);
    vt5meascpp=zeros(2,N);
    vt6meascpp=zeros(2,N);
    vt7meascpp=zeros(2,N);
    vt8meascpp=zeros(2,N);
    
    bsmeasoriginal=zeros(3,N);  % 先toa后dod最后aoa
    vt1measoriginal=zeros(3,N);
    vt2measoriginal=zeros(3,N);
    vt3measoriginal=zeros(3,N);
    vt4measoriginal=zeros(3,N);
    vt5measoriginal=zeros(3,N);
    vt6measoriginal=zeros(3,N);
    vt7measoriginal=zeros(3,N);
    vt8measoriginal=zeros(3,N);

    bsmeascpporiginal=zeros(2,N);  % 先toa后aoa
    vt1meascpporiginal=zeros(2,N);
    vt2meascpporiginal=zeros(2,N);
    vt3meascpporiginal=zeros(2,N);
    vt4meascpporiginal=zeros(2,N);
    vt5meascpporiginal=zeros(2,N);
    vt6meascpporiginal=zeros(2,N);
    vt7meascpporiginal=zeros(2,N);
    vt8meascpporiginal=zeros(2,N);

    gtPosetrans=gtPose';
    for i=1:N
        alpha=GetAngle(x(3:4,i),[0;0]);
        bsmeasoriginal(:,i)=[norm(x(1:2,i)-bs)+bias;getAOD([0;0;0],x(:,i),[0;0;0],'LOS');getAOA([0;0;0],x(:,i),[0;0;0],'LOS',alpha)'];
        vt1measoriginal(:,i)=[norm(x(1:2,i)-vt1(1:2))+vt1(3)+bias;getAOD([0;0;0],x(:,i),vt1,'VA');getAOA([0;0;0],x(:,i),vt1,'VA',alpha)'];
        vt2measoriginal(:,i)=[norm(x(1:2,i)-vt2(1:2))+vt2(3)+bias;getAOD([0;0;0],x(:,i),vt2,'VA')';getAOA([0;0;0],x(:,i),vt2,'VA',alpha)'];
        vt3measoriginal(:,i)=[norm(x(1:2,i)-vt3(1:2))+vt3(3)+bias;getAOD([0;0;0],x(:,i),vt3,'SP')';getAOA([0;0;0],x(:,i),vt3,'SP',alpha)'];
        vt4measoriginal(:,i)=[norm(x(1:2,i)-vt4(1:2))+vt4(3)+bias;getAOD([0;0;0],x(:,i),vt4,'VA')';getAOA([0;0;0],x(:,i),vt4,'VA',alpha)'];
        vt5measoriginal(:,i)=[norm(x(1:2,i)-vt5(1:2))+vt5(3)+bias;getAOD([0;0;0],x(:,i),vt3,'SP');getAOA([0;0;0],x(:,i),vt5,'SP',alpha)'];
        vt6measoriginal(:,i)=[norm(x(1:2,i)-vt6(1:2))+vt6(3)+bias;getAOD([0;0;0],x(:,i),vt3,'SP')';getAOA([0;0;0],x(:,i),vt6,'SP',alpha)'];
        vt7measoriginal(:,i)=[norm(x(1:2,i)-vt7(1:2))+vt7(3)+bias;getAOD([0;0;0],x(:,i),vt5,'SP')';getAOA([0;0;0],x(:,i),vt7,'SP',alpha)'];
        vt8measoriginal(:,i)=[norm(x(1:2,i)-vt8(1:2))+vt8(3)+bias;getAOD([0;0;0],x(:,i),vt6,'SP')';getAOA([0;0;0],x(:,i),vt8,'SP',alpha)'];
        
        bsmeascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-bs)+bias;GetAngle([0,0],gtPosetrans(2:3,i))];
        vt1meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt1(1:2))+vt1(3)+bias;GetAngle(vt1(1:2),gtPosetrans(2:3,i))];
        vt2meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt2(1:2))+vt2(3)+bias;GetAngle(vt2(1:2),gtPosetrans(2:3,i))];
        vt3meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt3(1:2))+vt3(3)+bias;GetAngle(vt3(1:2),gtPosetrans(2:3,i))];
        vt4meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt4(1:2))+vt4(3)+bias;GetAngle(vt4(1:2),gtPosetrans(2:3,i))];
        vt5meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt5(1:2))+vt5(3)+bias;GetAngle(vt5(1:2),gtPosetrans(2:3,i))];
        vt6meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt6(1:2))+vt6(3)+bias;GetAngle(vt6(1:2),gtPosetrans(2:3,i))];
        vt7meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt7(1:2))+vt7(3)+bias;GetAngle(vt7(1:2),gtPosetrans(2:3,i))];
        vt8meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt8(1:2))+vt8(3)+bias;GetAngle(vt8(1:2),gtPosetrans(2:3,i))];
    end
    
    for count=1:10
        % white noise addup
        bsnoise=[randn(1,N)*toalosrmse;randn(1,N)*aoalosrmse;randn(1,N)*aoalosrmse];
        vt1noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt2noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt3noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt4noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt5noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt6noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt7noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];
        vt8noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse;randn(1,N)*aoarmse];

        bsmeas=bsmeasoriginal+bsnoise;
        vt1meas=vt1measoriginal+vt1noise;
        vt2meas=vt2measoriginal+vt2noise;
        vt3meas=vt3measoriginal+vt3noise;
        vt4meas=vt4measoriginal+vt4noise;
        vt5meas=vt5measoriginal+vt5noise;
        vt6meas=vt6measoriginal+vt6noise;
        vt7meas=vt7measoriginal+vt7noise;
        vt8meas=vt8measoriginal+vt8noise;
        
        bsmeascpp=bsmeascpporiginal+bsnoise([1,3],:);
        vt1meascpp=vt1meascpporiginal+vt1noise([1,3],:);
        vt2meascpp=vt2meascpporiginal+vt2noise([1,3],:);
        vt3meascpp=vt3meascpporiginal+vt3noise([1,3],:);
        vt4meascpp=vt4meascpporiginal+vt4noise([1,3],:);
        vt5meascpp=vt5meascpporiginal+vt5noise([1,3],:);
        vt6meascpp=vt6meascpporiginal+vt6noise([1,3],:);
        vt7meascpp=vt7meascpporiginal+vt7noise([1,3],:);
        vt8meascpp=vt8meascpporiginal+vt8noise([1,3],:);

        bsmeascpp(2,:)=normalizeAngles(bsmeascpp(2,:));
        vt1meascpp(2,:)=normalizeAngles(vt1meascpp(2,:));
        vt2meascpp(2,:)=normalizeAngles(vt2meascpp(2,:));
        vt3meascpp(2,:)=normalizeAngles(vt3meascpp(2,:));
        vt4meascpp(2,:)=normalizeAngles(vt4meascpp(2,:));
        vt5meascpp(2,:)=normalizeAngles(vt5meascpp(2,:));
        vt6meascpp(2,:)=normalizeAngles(vt6meascpp(2,:));
        vt7meascpp(2,:)=normalizeAngles(vt7meascpp(2,:));
        vt8meascpp(2,:)=normalizeAngles(vt8meascpp(2,:));

        bsmeascpp(1,:)=normalizeRanges(bsmeascpp(1,:));
        vt1meascpp(1,:)=normalizeRanges(vt1meascpp(1,:));
        vt2meascpp(1,:)=normalizeRanges(vt2meascpp(1,:));
        vt3meascpp(1,:)=normalizeRanges(vt3meascpp(1,:));
        vt4meascpp(1,:)=normalizeRanges(vt4meascpp(1,:));
        vt5meascpp(1,:)=normalizeRanges(vt5meascpp(1,:));
        vt6meascpp(1,:)=normalizeRanges(vt6meascpp(1,:));
        vt7meascpp(1,:)=normalizeRanges(vt7meascpp(1,:));
        vt8meascpp(1,:)=normalizeRanges(vt8meascpp(1,:));

        h1=figure();
        hold on;
        plot(vt1meascpp(1,:));
        plot(vt2meascpp(1,:));
        plot(vt3meascpp(1,:));
        plot(vt4meascpp(1,:));
        plot(vt5meascpp(1,:));
        plot(vt6meascpp(1,:));
        plot(vt7meascpp(1,:));
        plot(vt8meascpp(1,:));

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
        para.N_VA = 3; % # number of VAs
        para.N_SP = 5;  % # number of SPs
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
            Channel.Values(:,6,i)=[vt5meas(1,i);vt5meas(2,i);0;vt5meas(3,i);0];
            Channel.Values(:,7,i)=[vt6meas(1,i);vt6meas(2,i);0;vt6meas(3,i);0];
            Channel.Values(:,8,i)=[vt7meas(1,i);vt7meas(2,i);0;vt7meas(3,i);0];
            Channel.Values(:,9,i)=[vt8meas(1,i);vt8meas(2,i);0;vt8meas(3,i);0];
            countall=-1;
            if i<=para.TIME/5
                v.Time(i).measurement(:,1)=Channel.Values(:,1,i);
                v.Time(i).measurement(:,2)=Channel.Values(:,2,i);
                v.Time(i).measurement(:,3)=Channel.Values(:,3,i);
                v.Time(i).measurement(:,4)=Channel.Values(:,4,i);
                v.Time(i).measurement(:,5)=Channel.Values(:,5,i);
                v.Time(i).measurement(:,6)=Channel.Values(:,6,i);
                v.Time(i).measurement(:,7)=Channel.Values(:,7,i);
                v.Time(i).measurement(:,8)=Channel.Values(:,8,i);
                v.Time(i).measurement(:,9)=Channel.Values(:,9,i);
                countall=9;
            else
                Channel.Visible.LOS(1,i)=0;
                Channel.Visible.TOT(1,i)=0;
                v.Time(i).measurement(:,1)=Channel.Values(:,2,i);
                v.Time(i).measurement(:,2)=Channel.Values(:,3,i);
                v.Time(i).measurement(:,3)=Channel.Values(:,4,i);
                v.Time(i).measurement(:,4)=Channel.Values(:,5,i);
                v.Time(i).measurement(:,5)=Channel.Values(:,6,i);
                v.Time(i).measurement(:,6)=Channel.Values(:,7,i);
                v.Time(i).measurement(:,7)=Channel.Values(:,8,i);
                v.Time(i).measurement(:,8)=Channel.Values(:,9,i);
                countall=8;
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
        VA(2).pos=[40; 0; 0];
        VA(3).pos=[40;10; 0];
        SP(1).pos=[5; -5; 0];
        SP(2).pos=[5; 15; 0];
        SP(3).pos=[15; -5; 0];
        SP(4).pos=[5; -5; 0];
        SP(5).pos=[5; -5; 0];

        if para.save ==1
            mkdir(['FinalData/Meas_',num2str(montecarlo)]);
            Sp = ['FinalData/Meas_',num2str(montecarlo),'/measurement_',num2str(count),'_' num2str(para.TIME)];
            filepath = Sp;
            save(filepath, 'state', 'para', 'BS', 'VA', 'SP', 'Channel', 'v', 'measurementCovariance')
            
            vtmeasurements=cell({vt1meascpp,vt2meascpp,vt3meascpp,vt4meascpp,vt5meascpp,vt6meascpp,vt7meascpp,vt8meascpp});
            % saving
            meas_complete=cell(N,1);

            N_clutter = poissrnd(lambda,1);

            for i=1:N
                meas_complete{i}=cell(0,1);
                mp = struct(OriginPosition = [gtPose(i,2:3), bias]);
                for j=1:8
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

            target=[vt1,vt2,vt3,vt4,vt5,vt6,vt7,vt8];
            save(['FinalData/Meas_',num2str(montecarlo),'/meas_',num2str(count),'_phd.mat'],'meas_complete','target','gtPose','bsmeascpp');
            print (h1, '-dpng', ['FinalData/Meas_',num2str(montecarlo),'/meas_',num2str(count),'.png']);
            close(h1);
        end
        print (h, '-dpng', ['FinalData/Meas_',num2str(montecarlo),'/gtTraj.png']);
        save( ['FinalData/Meas_' num2str(montecarlo) '/gtPose.mat'],'gtPose');
        clear parameter;
        clear Channel;
        clear v;
        
    end
end

