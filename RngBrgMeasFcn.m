function [rngbrg] = RngBrgMeasFcn(x)
%RNGBRG 此处显示有关此函数的摘要
%   此处显示详细说明
[sizer,sizec]=size(x);
rngbrg=zeros(2,sizec);
rngbrg(1,:)=sqrt(x(1,:).^2+x(2,:).^2);
for i=1:sizec
    rngbrg(2,i)=standardatan(x(1,i),x(2,i));
end
end

