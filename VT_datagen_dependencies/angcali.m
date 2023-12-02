function [out] = angcali(in)
%UNTITLED 此处显示有关此函数的摘要
%   -pi~pi
out=in;
while(out>pi)
    out=out-2*pi;
end
while(out<-pi)
    out=out+2*pi;
end
end

