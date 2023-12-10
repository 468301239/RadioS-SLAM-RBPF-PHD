function [weightfactor] = multifeaturereweight(measures,filter,x,maxnumber)
P_d=0.9;
[~, sortedIndices] = sort(filter.phd.Weights, 'descend');
Evaln=min(length(extractState(filter.phd,0.007)),maxnumber);
EvalStates=filter.phd.States(:,sortedIndices(1:Evaln));
EvalStateCovs=filter.phd.StateCovariances(:,:,sortedIndices(1:Evaln));
if Evaln<=0
    weightfactor=1;
    return;
end
vbefore=1;
vafter=1;
for i=1:Evaln
    vbefore=vbefore*calculatePHD(filter.phd_old,EvalStates(:,i));
    vafter=vafter*calculatePHD(filter.phd,EvalStates(:,i));
end
md2threshold=50;
nZ=length(measures);
likelihoodmatrix=zeros(Evaln,nZ);
for i=1:Evaln
    for j=1:nZ
        z=measures{j}.Measurement;
        mp = struct(OriginPosition = x);
        zexp=RngBrgMeasFcnVT(EvalStates(:,i),mp);
        zJac=RngBrgMeasFcnVTjac(EvalStates(:,i),mp);
        stateCov=EvalStateCovs(:,:,i);
        zCov=zJac*stateCov*zJac'+diag([filter.config.Rngnoise^2,filter.config.Brgnoise^2]);
        error = z-zexp;
        error(2)=normalizeAngles(error(2));
        md2=error'/zCov*error;
        if md2<md2threshold
            likelihoodmatrix(i,j)=1/sqrt(det(2*pi*zCov))*exp(-1/2*md2);
        end
    end
end
assignmentlist=zeros(nZ,1);
assignmentlikelihood=zeros(nZ,1);
measweight=1;
for i=1:nZ
    [assignmentlikelihood(i),assignmentlist(i)]=max(likelihoodmatrix(:,i));
    if assignmentlikelihood(i)>0
        measweight=measweight*assignmentlikelihood(i)*P_d;
    else
        measweight=measweight*(1-P_d)*filter.config.Kc;
    end
end
weightfactor=measweight*exp(filter.sumweight-filter.sumweight_old)*vbefore/vafter;
end