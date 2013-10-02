function PlotExampleSpectrum(mod_struct,M,osr,f0,quadrature)
%PlotExampleSpectrum(ntf|mod_struct,M=1,osr=64,f0=0,quadrature=0)
% ntf|mod_struct is either the NTF or a struct containing 
% ntf, M=1, osr=64, f0=0, quadrature=0

% Handle the input arguments
parameters = {'mod_struct' 'M' 'osr' 'f0' 'quadrature'};
defaults = { [] 1 64 0 0 };
for arg_ii=1:length(defaults)
    parameter = parameters{arg_ii};
    if arg_ii>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_ii};'])
    end
end
if nargin==1 & isstruct(mod_struct)
    flattenStruct(mod_struct);
else
    ntf = mod_struct;
end
[f1 f2] = ds_f1f2(osr,f0,quadrature);

delta = 2;

% Plot an example spectrum
Amp = undbv(-3);  % Test tone amplitude, relative to full-scale.
f = 0.3;          % Test tone frequency offset from f0, relative to bw. 
                  % (Will be adjusted to be an fft bin)
N = 2^12;
f1_bin = round(f1*N);
f2_bin = round(f2*N);
fin = round(((1-f)/2*f1 + (f+1)/2*f2)*N);
if ~quadrature
    u = Amp*M*cos((2*pi/N)*fin*[0:N-1]);
    v = simulateDSM(u,ntf,M+1);
else
    u = Amp*M*exp((2i*pi/N)*fin*[0:N-1]);
    v = simulateQDSM(u,ntf,M+1);
end
window = ds_hann(N);
NBW = 1.5/N;
spec0 = fft(v.*window)/(M*N/4);
if ~quadrature
    freq = linspace(0,0.5,N/2+1);
	plot(freq,dbv(spec0(1:N/2+1)),'c','Linewidth',1);
	hold on
	spec_smoothed = circ_smooth(abs(spec0).^2, 16);
	plot(freq, dbp(spec_smoothed(1:N/2+1)), 'b', 'Linewidth', 3);
	Snn = abs(evalTF(ntf,exp(2i*pi*freq))).^2 * 2/12 * (delta/M)^2;
	plot(freq, dbp(Snn*NBW), 'm', 'Linewidth', 1);
    snr = calculateSNR(spec0(f1_bin+1:f2_bin+1),fin-f1_bin);
    msg = sprintf('SQNR = %.1fdB\n @ A=%.1fdBFS & osr=%.0f\n', ...
      snr, dbv(spec0(fin+1)), osr );
    if f0<0.25
        text(f0+1/osr, -15, msg,'Hor','Left');
    else
        text(f0-1/osr, -15, msg,'Hor','Right');
    end
    text(0.5,-135,sprintf('NBW=%.1e ',NBW),'hor','right');
    figureMagic([0 0.5],1/16,4, [-140 0],10,2);
else
    spec0 = fftshift(spec0/2);
    freq = linspace(-0.5,0.5,N+1); freq(end)=[];
	plot(freq,dbv(spec0),'c','Linewidth',1);
	hold on
	spec_smoothed = circ_smooth(abs(spec0).^2, 16);
	plot(freq, dbp(spec_smoothed), 'b', 'Linewidth', 3);
	Snn = abs(evalTF(ntf,exp(2i*pi*freq))).^2 * 2/12 * (delta/M)^2;
	plot(freq, dbp(Snn*NBW), 'm', 'Linewidth', 1);
    snr = calculateSNR(spec0(N/2+1+[f1_bin:f2_bin]),fin-f1_bin);
    msg = sprintf('SQNR = %.1fdB\n @ A=%.1fdBFS & osr=%.0f\n', ...
      snr, dbv(spec0(N/2+fin+1)), osr );
    if f0>=0
        text(f0-0.05, -15, msg,'Hor','Right');
    else
        text(f0+0.05, -15, msg,'Hor','Left');
    end
    text(-0.5,-135,sprintf(' NBW=%.1e',NBW),'hor','left');
    figureMagic([-0.5 0.5],0.125,2, [-140 0],10,2);
end
xlabel('frequency');
% ylabel('PSD (dBFS/NBW)');
% printmif(fullfile('MIF','spec'), [5 2], 'Helvetica10')
