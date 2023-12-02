function [Jmx] = RngBrgMeasFcnjac(x)
%RNGBRGMEASFCNJAC 此处显示有关此函数的摘要
%   此处显示详细说明
dis=sqrt(x(1)^2+x(2)^2);
Jmx=[
x(1)/dis, x(2)/dis;
-x(2)/dis^2, x(1)/dis^2
];
end

