function [out] = GetAngle(x1,x2)
%AngleGet return [-pi,pi] x1相对于x2
x=x1-x2;
if(x(1)==0)
    if(x(2)>=0)
        out=pi/2;
    else
        out=-pi/2;
    end
else
    out=atan(x(2)/x(1));
    if(x(1)<0)
        out=out+pi;
    end
    if(out<-pi)
        out=out+2*pi;
    end
    if(out>pi)
        out=out-2*pi;
    end
end
end