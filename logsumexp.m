function s = logsumexp(x)
% This function calculates logsumexp(x). 
% logsumexp(x) = log(x(1) + x(2) + x(3) ..) given log(x(1)), log(x(2)) .. 
xmax = max(x);
xdiff = x - xmax;
s = xmax + log(sum(exp(xdiff)));
end