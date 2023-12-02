function [rngbrg] = RngBrgMeasFcnVT(x,parameters)
%RNGBRG 此处显示有关此函数的摘要
%   此处显示详细说明
temp=parameters.OriginPosition;
orix=temp(1);
oriy=temp(2);
orib=temp(3);
[sizer,sizec]=size(x);
rngbrg=zeros(2,sizec);
rngbrg(1,:)=sqrt((x(1,:)-orix).^2+(x(2,:)-oriy).^2)+orib+x(3,:);
for i=1:sizec
    rngbrg(2,i)=standardatan((x(1,i)-orix),(x(2,i)-oriy));
end
end

