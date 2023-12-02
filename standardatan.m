function [theta] = standardatan(x,y)
%STANDARDATAN 此处显示有关此函数的摘要
%   theta \in [-pi,pi]
if x>0
    theta=atan(y/x);
elseif x<0
    if y>0
        theta=atan(y/x)+pi;
    else
        theta=atan(y/x)-pi;
    end
else
    theta=(double(y>0)-1/2)*pi;
end
end