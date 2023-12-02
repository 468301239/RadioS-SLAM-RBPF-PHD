clear all; close all; clc;
f=12.5;   % frequency hz
dt=1/f; % period s
T=30;   % exp time s
N=floor(T/dt);

bs=[0;0];
vt1=[0;10;0];
vt2=[10;-5;sqrt(125)];
vt3=[10;15;sqrt(125)];
vt4=[10;-5;sqrt(325)];

ax=0.5; ay=ax;    % acceleration for each direction m/s^2
pos_init=[0;0]; % initial position m
v_init=[1;0];   % initial velocity m/s

toarmse=0.3; % m
toalosrmse=0.1; % m
aoarmse=5*pi/180;   % rad
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
for montecarlo=1:1
    flag=0;
    while (flag~=1)
        a=[normrnd(0,1,1,N)*ax;normrnd(0,1,1,N)*ay];
        for i=2:N
            x(:,i)=A*x(:,i-1)+B*a(:,i);
            if(x(2,i)>5)
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

    bsmeascpp=zeros(2,N);  % 先toa后dod最后aoa
    vt1meascpp=zeros(2,N);
    vt2meascpp=zeros(2,N);
    vt3meascpp=zeros(2,N);
    vt4meascpp=zeros(2,N);

    bsmeascpporiginal=zeros(2,N);  % 先toa后dod最后aoa
    vt1meascpporiginal=zeros(2,N);
    vt2meascpporiginal=zeros(2,N);
    vt3meascpporiginal=zeros(2,N);
    vt4meascpporiginal=zeros(2,N);

    gtPosetrans=gtPose';
    for i=1:N
        bsmeascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-bs)+bias;GetAngle([0,0],gtPosetrans(2:3,i))];
        vt1meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt1(1:2))+vt1(3)+bias;GetAngle(vt1(1:2),gtPosetrans(2:3,i))];
        vt2meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt2(1:2))+vt2(3)+bias;GetAngle(vt2(1:2),gtPosetrans(2:3,i))];
        vt3meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt3(1:2))+vt3(3)+bias;GetAngle(vt3(1:2),gtPosetrans(2:3,i))];
        vt4meascpporiginal(:,i)=[norm(gtPosetrans(2:3,i)-vt4(1:2))+vt4(3)+bias;GetAngle(vt4(1:2),gtPosetrans(2:3,i))];
    end

    % white noise addup
    bsnoise=[randn(1,N)*toalosrmse;randn(1,N)*aoalosrmse];
    vt1noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse];
    vt2noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse];
    vt3noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse];
    vt4noise=[randn(1,N)*toarmse;randn(1,N)*aoarmse];

    bsmeascpp=bsmeascpporiginal+bsnoise;
    vt1meascpp=vt1meascpporiginal+vt1noise;
    vt2meascpp=vt2meascpporiginal+vt2noise;
    vt3meascpp=vt3meascpporiginal+vt3noise;
    vt4meascpp=vt4meascpporiginal+vt4noise;
    
    bsmeascpp(2,:)=normalizeAngles(bsmeascpp(2,:));
    vt1meascpp(2,:)=normalizeAngles(vt1meascpp(2,:));
    vt2meascpp(2,:)=normalizeAngles(vt2meascpp(2,:));
    vt3meascpp(2,:)=normalizeAngles(vt3meascpp(2,:));
    vt4meascpp(2,:)=normalizeAngles(vt4meascpp(2,:));

    vtmeasurements=cell({vt1meascpp,vt2meascpp,vt3meascpp,vt4meascpp});
    % saving
    meas_complete=cell(N,1);
        
    N_clutter = poissrnd(lambda,1);
    
    for i=1:N
        meas_complete{i}=cell(0,1);
        mp = struct(OriginPosition = [gtPose(i,2:3), bias]);
        for j=1:4
            if rand()<P_D
                tempcell=cell(1,1);
                tempcell{1,1}=objectDetection((i-1)*dt,vtmeasurements{j}(:,i),'MeasurementNoise',R,MeasurementParameters=mp);
                meas_complete{i}=[meas_complete{i};tempcell];
            end
        end
        
        for j=1:N_clutter
            r_clutter=50*rand();
            b_clutter=2*pi*rand()-pi;
            tempcell=cell(1,1);
            tempcell{1,1}=objectDetection((i-1)*dt,[r_clutter;b_clutter],'MeasurementNoise',R,MeasurementParameters=mp);
            meas_complete{i}=[meas_complete{i};tempcell];
        end
    end
    
    target=[vt1,vt2,vt3,vt4];
    print (h, '-dpng', 'dataVT/test.png');
    save('dataVT/meas.mat','meas_complete','target','gtPose','bsmeascpp');
end