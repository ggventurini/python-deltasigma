function [sv,sx] = simulateBiDWA(v,M,sx)
% [sv,sx] = simulateBiDWA(v,M,sx=0)    Simulate the Bi-directional DWA algorithm
% I. Fujimori, L. Longo, A. Hairapetian, K. Seiyama, S. Kosic, J.  Cao and S.L. Chan, 
% "A 90-dB SNR 2.5-MHz output-rate ADC using cascaded multibit delta-sigma modulation at 8x oversampling ratio,"
% IEEE Journal of Solid-State Circuits, vol.35, no.12, pp.1820-1828, Dec. 2000.

% Argument checking and default-setting 
ArgumentsAndDefaults = {
 'v'   NaN
 'M'   16
 'sx'  [0 0  0]
  };
for i = 1:size(ArgumentsAndDefaults,1)
    parameter = ArgumentsAndDefaults{i,1};
    if i>nargin || eval(['isempty(' parameter ') ']) || ...
      ( eval(['isnumeric(' parameter ') '])  &&  ...
        eval(['length(' parameter ') <= 1']) && ...
        eval(['isnan(' parameter ')']))
        if isnan(ArgumentsAndDefaults{i,2})
            error('%s: Argument %d (%s) is required.',mfilename, i, parameter )
        else
            eval([parameter '= ArgumentsAndDefaults{i,2};'])
        end
    end
end

N = length(v);
p0 = sx(1);
p1 = sx(2);
dir = sx(3);

v = (v+M)/2; % Translate -M:2:M to 0:M
sv = zeros(M,N);
for n=1:N
    if dir==0
        pp = p0+v(n);
        i = mod(p0:pp-1,M);
        sv(i+1,n) = 1;
        p0 = mod(pp,M);
        dir = 1;
    else
        pp = p1-v(n);
        i = mod(pp:p1-1,M);
        sv(i+1,n) = 1;
        p1 = mod(pp,M);
        dir = 0;
    end
end

sv = 2*sv - 1; % Translate [0 1] to [-1 1]
sx = [p0 p1 dir];


