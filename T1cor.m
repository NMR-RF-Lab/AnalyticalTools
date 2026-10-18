function [Ac, Acstd] = T1cor(T1, T1std, TR, alpha, A, Astd)
% T1COR Computes the T1-corrected amplitude and its propagated uncertainty.
%
% INPUTS:
%   A         - Measured peak amplitude
%   sigma_A   - Standard deviation of peak amplitude
%   T1        - T1 relaxation time (same units as TR)
%   sigma_T1  - Standard deviation of T1
%   TR        - Repetition time (same units as T1)
%   alpha     - Flip angle in radians
%
% OUTPUTS:
%   Ac        - T1-corrected amplitude
%   sigma_Ac  - Propagated uncertainty (standard deviation) of corrected amplitude

% --- Precompute common terms ---
E        = exp(-TR./T1);
cos_a    = cosd(alpha);

% --- Corrected amplitude ---
Ac = A.*(1 - cos_a.*E)./(1 - E);

% --- Partial derivative: dAc/dA ---
dAc_dA = (1 - cos_a.*E)./(1 - E);

% --- Partial derivative: dAc/dT1 ---
dAc_dT1 = A.*(1 - cos_a).*E.*TR./(T1.^2.*(1 - E).^2);

% --- Propagated uncertainty ---
if nargin == 6
    Acstd = sqrt((dAc_dA.*Astd).^2 + (dAc_dT1.*T1std).^2);
end
end