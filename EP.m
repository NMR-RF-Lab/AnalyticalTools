function ERR = EP(type, err, cv, data)
% Error propagation
% Inputs:
%   - err: uncertainty, standard deviation or CRLB (in std units)
%   - cv: covariance
switch lower(type)
    case 'addition'
        ERR = sqrt(sum(err.^2,[],"omitmissing") + 2*sum(cv,[],"omitmissing"));
    case 'subtraction'
        ERR = sqrt(sum(err.^2,[],"omitmissing") - 2*sum(cv,[],"omitmissing"));
    case 'multiplication'
        ERR = abs(prod(data,[],"omitmissing"))*...
            sqrt(sum((err./data).^2,[],"omitmissing") + 2*sum(cv./(data(2:end).*data(1:end-1)),[],"omitmissing"));
    case 'division'
        ERR = abs(prod(data(1:end-1),[],"omitmissing")/data(end))*...
            sqrt(sum((err./data).^2,[],"omitmissing") - 2*sum(cv./(data(2:end).*data(1:end-1)),[],"omitmissing"));
end
end