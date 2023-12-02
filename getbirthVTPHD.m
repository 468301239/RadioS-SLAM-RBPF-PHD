function [birthPHD] = getbirthVTPHD(reservedMeasurements,wb,R)
%GETBIRTHPHD 此处显示有关此函数的摘要
%   此处显示详细说明
    gamma=0.7;
    zeta=0.1;
    iota=0.5;
    xi=0.3;
    mb = zeros(3,length(reservedMeasurements));
    Pb = repmat(zeros(3),1,1,length(reservedMeasurements));
    for i=1:length(reservedMeasurements)
        Rng=reservedMeasurements{i}.Measurement(1);
        Brg=reservedMeasurements{i}.Measurement(2);
        temp=reservedMeasurements{i}.MeasurementParameters.OriginPosition;
        orix=temp(1);
        oriy=temp(2);
        orib=temp(3);
        mb(1,i)=orix+(Rng-orib)*cos(Brg)*gamma;
        mb(2,i)=oriy+(Rng-orib)*sin(Brg)*gamma;
        mb(3,i)=(Rng-orib)*(1-gamma);
        T=[cos(Brg)/sqrt(2),-sin(Brg),cos(Brg)/sqrt(2);
            sin(Brg)/sqrt(2),cos(Brg),sin(Brg)/sqrt(2);
            -1/sqrt(2),0,1/sqrt(2)
        ];
        S=diag([zeta*(Rng-orib)^2,iota*(Rng-orib)^2*R(2,2),xi*R(1,1)]);
        Pb(:,:,i)=T*S*T';
    end
    wblist=wb*ones(1,length(reservedMeasurements));
    birthPHD = gmphd(mb,Pb,'Weights',wblist,'MaxNumComponents',20);
end

