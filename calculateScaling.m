function Lij = calculateScaling(logqij, w, Pd, Kc)
% Our objective is to calculate
% Lij = Pd*qij/(Kc + sumOveri(Pd*qij*wi));
logPdqij = log(Pd) + logqij;

logPdqijwi = log(Pd) + logqij + log(w(:));

% Allocate sum of logPdqijwi over i.
logSumOverIPdqijwi = zeros(1,size(logqij,2));
for j = 1:size(logqij,2)
    logSumOverIPdqijwi(j) = logsumexp(logPdqijwi(:,j));
end

% Compute log(Kc + sumOveri(wi*Pd*qij));
den = zeros(1,size(logqij,2));
logKc = log(Kc);
for j = 1:size(logqij,2)
    den(j) = logsumexp([logKc;logSumOverIPdqijwi(j)]);
end

logLij = logPdqij - den;
Lij = exp(logLij);

end