function [Ac, Acstd] = STEAMrelaxcor(T1, T1std, T2, T2std, TR, TE, TM, A, Astd)
% STEAMRELAXCOR computes relaxation-corrected amplitude for a STEAM
% sequence. This assumes a perfect 90 degree excitation. This also assumes
% independent uncertainties between A, T1, and T2.
%
% INPUTS:
%   T1      - T1 relaxation time, same units as TR and TM
%   T1std   - Standard deviation of T1
%   T2      - T2 relaxation time, same units as TE
%   T2std   - Standard deviation of T2
%   TR      - Repetition time, same units as T1
%   TE      - Echo time, same units as T2
%   TM      - STEAM mixing time, same units as T1
%   A       - Measured peak amplitude or area
%   Astd    - Standard deviation of measured peak amplitude or area
%
% OUTPUTS:
%   Ac      - Relaxation-corrected amplitude
%   Acstd   - Propagated uncertainty

% --- Exponential recovery/decay terms ---
E = exp(-TR./T1);

% --- Combined correction factor ---
F = exp(TE./T2).*exp(TM./T1)./(1 - E);

% --- Corrected amplitude ---
Ac = A.*F;

% --- Partial derivative with respect to A ---
dAc_dA = Ac./A;

% --- Partial derivative with respect to T1 ---
dAc_dT1 = Ac.*(-TM./(T1.^2) + TR.*E./((T1.^2).*(1 - E)));

% --- Partial derivative with respect to T2 ---
dAc_dT2 = Ac.*-TE./(T2.^2);

% --- Propagated uncertainty ---
if nargin == 9
    Acstd = sqrt((dAc_dA.*Astd).^2 + (dAc_dT1.*T1std).^2 + (dAc_dT2.*T2std).^2);
end
end