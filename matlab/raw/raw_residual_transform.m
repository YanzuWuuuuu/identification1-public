function transformed=raw_residual_transform(signal,ts,band)
% Fixed residual weighting, not nonparametric system estimation.
if isempty(band), transformed=signal; return; end
N=size(signal,1); window=.5-.5*cos(2*pi*(0:N-1)'/(N-1));
frequencies=(0:N-1)'/(N*ts);
keep=frequencies>=band(1) & frequencies<=band(2);
spectrum=fft(signal.*window,[],1);
transformed=[real(spectrum(keep,:));imag(spectrum(keep,:))];
end
