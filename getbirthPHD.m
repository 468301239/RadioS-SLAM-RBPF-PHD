function [birthPHD] = getbirthPHD(reservedMeasurements,wb,R)
%GETBIRTHPHD 此处显示有关此函数的摘要
%   此处显示详细说明
    mb = zeros(2,length(reservedMeasurements));
    Pb = repmat(zeros(2),1,1,length(reservedMeasurements));
    for i=1:length(reservedMeasurements)
        Rng=reservedMeasurements{i}.Measurement(1);
        Brg=reservedMeasurements{i}.Measurement(2);
        mb(1,i)=Rng*cos(Brg);
        mb(2,i)=Rng*sin(Brg);
        H=[cos(Brg),-Rng*sin(Brg);sin(Brg),Rng*cos(Brg)];
        Pb(:,:,i)=H*R*H';
    end
    wblist=wb*ones(1,length(reservedMeasurements));
    birthPHD = gmphd(mb,Pb,'Weights',wblist,'MaxNumComponents',20);
end

