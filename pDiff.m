% PDIFF calculates the percent difference between a new value (A) and an
% old value (B). If B == 0, then the percent difference is set to 0.
function pd = pDiff(A,B)
if B == 0
    pd = 0;
else
    pd = abs(A-B)/B;
end
end