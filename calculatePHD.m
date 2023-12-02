function D = calculatePHD(phdObject, x)
D = 0;
m = phdObject.States;
P = phdObject.StateCovariances;
w = phdObject.Weights;
for i = 1:phdObject.NumComponents
    e = m(:,i) - x;
    gaussianPdf = 1/sqrt(det(2*pi*P(:,:,i)))*exp(-1/2*e'/P(:,:,i)*e);
    D = D + w(i)*gaussianPdf;
end
end
