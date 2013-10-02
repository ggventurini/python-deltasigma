function w = ds_hann(n)
% function w = ds_hann(n)
% A Hann window of length n. Does not smear tones located exactly in a bin.

% Note: This function was formerly just "hann." Re-naming
% was necessary to avoid a conflict with Mathworks's function
% of the same name.
w = .5*(1 - cos(2*pi*(0:n-1)/n) );
