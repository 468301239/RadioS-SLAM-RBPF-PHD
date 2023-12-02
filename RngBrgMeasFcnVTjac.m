function [Jmx] = RngBrgMeasFcnVTjac(x,parameters)
%RNGBRGMEASFCNJAC 此处显示有关此函数的摘要
%   此处显示详细说明
temp=parameters.OriginPosition;
orix=temp(1);
oriy=temp(2);
orib=temp(3);
dis=sqrt((x(1)-orix)^2+(x(2)-oriy)^2);
Jmx=[
(x(1)-orix)/dis, (x(2)-oriy)/dis,   1;
-(x(2)-oriy)/dis^2, (x(1)-orix)/dis^2,  0
];
end

