clear;
clc;
lim=40;
target=[-40;0];
target=[4,4,-4,-4,10,10,-10,-10,16,16,-16,-16,40,40,-40,-40,-30;4,-10,16,-40,4,-10,16,-40,4,-10,16,-40,4,-10,16,-40,0];
Rngnoise=0.3;
Brgnoise=5*pi/180;
R=diag([Rngnoise^2,Brgnoise^2]);
Pd=0.95;

N=200;
meas_complete=cell(N,1);

for i=1:N
    meas_complete{i}=cell(0,1);
    for j=1:17
        if rand()<Pd
            w=diag([Rngnoise,Brgnoise])*randn(2,1);
            meas_temp=RngBrgMeasFcn(target(:,j))+w;
            meas_temp(2,:)=normalizeAngles(meas_temp(2,:));
            tempcell=cell(1,1);
            tempcell{1,1}=objectDetection(i,meas_temp,'MeasurementNoise',R);
            meas_complete{i}=[meas_complete{i};tempcell];
        end
    end
end

h=figure();
scatter(0,0);
hold on;
for i=1:N
    for j=1:length(meas_complete{i})
        scatter(meas_complete{i}{j}.Measurement(1)*cos(meas_complete{i}{j}.Measurement(2)),meas_complete{i}{j}.Measurement(1)*sin(meas_complete{i}{j}.Measurement(2)));
    end
end

hold on;
box on;

scatter(target(1,:),target(2,:),'Marker','s','MarkerFaceColor','b','MarkerEdgeColor','none')

print (h, '-dpng', 'data/test.png');
save('data/meas.mat','meas_complete','target');