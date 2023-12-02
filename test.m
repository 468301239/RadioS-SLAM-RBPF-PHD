timeindex=20;
gridnum=50;
origpos=meas_complete{timeindex}{1}.MeasurementParameters.OriginPosition;
x=linspace(origpos(1)-6,origpos(1)+6,gridnum);
y=linspace(origpos(2)-6,origpos(2)+6,gridnum);
[X, Y] = meshgrid(x, y);
weighting=zeros(gridnum,gridnum);
for i=1:gridnum
    for j=1:gridnum
        weighting(i,j)=multifeaturereweight(meas_complete{timeindex},myfilter,[X(i,j),Y(i,j),0.3],100);
    end
end

figure;
surf(X, Y, weighting);
title('Three-dimensional Plot');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');