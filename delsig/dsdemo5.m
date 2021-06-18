% Demonstrate the simulateMS function
echo off;
if ~exist('LiveDemo','var')
    LiveDemo = 0;
end

show_usage = 1;
osr = 25;
ntf = synthesizeNTF(6,osr,1,4);
M = 16;
N = 2^14;
sigma_d = 0.01;		% 1% mismatch
window = ds_hann(N)/(M*N/8);
windowM = repmat(window,M,1);

mtf1 = zpk(1,0,1,1);                    % First-order shaping
mtf2 = zpk([ 1 1 ], [ 0.3 0.3 ], 1, 1);	% Second-order shaping
mtf2b = synthesizeNTF(2,osr,1,2);	    % Second-order shaping
mtf4 = synthesizeNTF(4,osr*0.9,1,1.3);	% Fourth-order shaping
dw1 = [8 8 4 4 2 2 1 1 zeros(1,8)]';

cases = {
    'A'        'f'  'mtf' 'dither' 'dw'  'leg_i'        %case number
    undbv(-3)  0.01   []      0     []   'Thermometer'  %1
    undbv(-3)  0.01  mtf1     0     []   'Rotation'     %2
    undbv(-3)  0.01  mtf2     0     []   '2^{nd}-order' %3
    undbv(-30) 0.01  mtf1     0     []   'Rotation'     %4
    undbv(-30) 0.01  mtf1    0.5    []   'Rot + dither'    %5
    undbv(-3)  0.01  mtf2b    0     []   '2^{nd}-order with zero'    %6
    undbv(-3)  0.01  mtf4     0     []   '4^{th}-order' %7
   };
comparisons = { [1 2], [2 3], [4 5], [3 6], [6 7] };

%% Run simulations
sv = cell(size(cases,1)-1,1);
leg = sv;
Svv = sv;
Sdd = sv;
fin = zeros(size(cases,1)-1,1);
for i = 1:size(cases,1)-1
    for j = 1:size(cases,2) % Set the variables for each case
        eval([cases{1,j} '= cases{i+1,j};'])
    end
    fin(i) = round(f*N);
    inband = setdiff(2:1+ceil(0.5*N/osr),1+[0 1 fin(i)+[-1 0 1]]);
    w = (2*pi/N)*fin(i);
    u = M*A*sin(w*(0:N-1));
    v = simulateDSM(u,ntf,M+1);	% M unit elements requires an (M+1)-level quant.
    Svv{i} = abs(fft(v.*window)).^2;
    if isempty(mtf)
        sv{i} = ds_therm(v,M);
    else
        sv{i} = simulateMS(v,M,mtf,dither,dw);
    end
    Sdd{i} = sigma_d^2 * sum( abs(fft(sv{i}.*windowM,[],2)).^2 );
    mnp = sum(Sdd{i}(inband))/1.5;
    leg{i} = sprintf('%s (MNP= %.0f dBFS)',leg_i,dbp(mnp));
end

%% Plot results
 for comp_i=1:length(comparisons)
    clc
    case_nums = comparisons{comp_i};
    fprintf(1,'\t\tMismatch-Shaping Unit-Element DAC\n\n');
    fprintf(1,'Comparing %s vs. %s.\n', leg{case_nums(1)}, leg{case_nums(2)});
    nc = length(case_nums);
    if show_usage
        figure(1); clf
        set(gcf,'NumberTitle','off');
        set(gcf,'Name','Element Usage');
        T = 25;
        for i = 1:nc
            subplot(nc,1,i);
            ci = case_nums(i);
            plotUsage(sv{ci}(:,1:T));
        end
    end
    if LiveDemo
        set(1,'position',[9 204 330 525]);
        changeFig(18,.5,1);
        pause
    end

    figure(2); clf;
    set(gcf,'NumberTitle','off');
    set(gcf,'Name','Error Spectra');
    cols = {'b1','m1','r1'};
    cleg = cell(nc,1);
    for i = 1:nc
        ci = case_nums(i);
        plotSpectrum(sqrt(Sdd{ci}),fin(ci),cols{i},[],4);
        cleg{i} = leg{ci};
        hold on;
    end
    plotSpectrum(sqrt(Svv{ci}),fin(ci),'g',[],5);
    axis([1e-3 0.5 -140 -50]);
    plot([1e-3 0.5/osr],-140*[1 1],'-k','Linewidth',4);
    text(0.5,-140,sprintf('NBW=%.1e ',1.5/N),'Hor','right','Ver','bot');
    grid on;
    ylabel('Error PSD');
    xlabel('Normalized Frequency');
    title(sprintf('A = %.0fdBFS', dbp(Svv{ci}(fin(ci)))));
    legend(cleg,'Location','Northeast');
    if LiveDemo
        figure(1);
        set(1,'position',[9 427 200 300]);
        changeFig;
        figure(2);
        set(2,'position',[237 288 553 439]);
        changeFig(18,2,12);
        pause
        changeFig;
        legend;
        set(2,'position',[238 427 325 300]);
    elseif comp_i<length(comparisons)
        fprintf(1,'Paused.\n');
        pause
    end
end

if 0 % Quadrature example
    order = 4;
    osr = 32;
    M = 8;
    f0 = 1/16;
    quadrature = 1;
    [f1 f2] = ds_f1f2(osr,f0,quadrature);
    ntf = synthesizeQNTF(order,osr,f0,-50,-10);
    Amp = undbv(-3);  % Test tone amplitude, relative to full-scale.
    f = 0.3;          % Test tone frequency offset from f0, relative to bw. 
                      % (Will be adjusted to be an fft bin)
    N = 2^12;
    f1_bin = round(f1*N);
    f2_bin = round(f2*N);
    fin = round(((1-f)/2*f1 + (f+1)/2*f2)*N);

    u = Amp*M*exp((2i*pi/N)*fin*(0:N-1));
    v = simulateQDSM(u,ntf,M+1);
    sv0 = [1i*ds_therm(imag(v),M); 
           ds_therm(real(v),M)];
    mtf1 = zpk(exp(2i*pi*f0),0,1,1);	%First-order complex shaping
    sv1 = simulateQESL(v,mtf1,M);
    figure(1); clf
    subplot(121);
    plotUsage(sv0(:,1:20));
    subplot(122);
    plotUsage(sv1(:,1:20));
    %printmif(fullfile('MIF','QUsage'), [32 4], 'Helvetica8')    
    ue = 1 + 0.01*randn(2*M,1);
    dv0 = ue' * sv0;
    dv1 = ue' * sv1;
    window = ds_hann(N);
    Nsm = 32;
    NBW = 1.5/N * Nsm/2 / 1.5;  % Approximate correction for smoothing; amplitude correction is also approximate
    spec =  fft(v.*window)/(M*N/2);
    spec0 = fft(dv0.*window)/(M*N/2);
    spec1 = fft(dv1.*window)/(M*N/2);
    spec_sm = circ_smooth(abs(fftshift(spec)).^2,Nsm)* Nsm/3;
    spec0_sm = circ_smooth(abs(fftshift(spec0)).^2,Nsm)* Nsm/3;
    spec1_sm = circ_smooth(abs(fftshift(spec1)).^2,Nsm)* Nsm/3;
    freq = linspace(-0.5,0.5,N+1); freq(end)=[];
    figure(2); clf
	plot(freq,dbp(spec_sm),'k','Linewidth',1);
    hold on;
	plot(freq,dbp(spec0_sm),'r','Linewidth',1);
	plot(freq,dbp(spec1_sm),'b','Linewidth',1);
    snr = calculateSNR(spec(f1_bin+1:f2_bin+1),fin-f1_bin);
    snr0 = calculateSNR(spec0(f1_bin+1:f2_bin+1),fin-f1_bin);
    snr1 = calculateSNR(spec1(f1_bin+1:f2_bin+1),fin-f1_bin);
    msg = sprintf('Ideal: SNR = %.0fdB', snr);
    msg0 = sprintf('No MS: SNR = %.0fdB', snr0);
    msg1 = sprintf('1^{st}-order MS: SNR = %.0fdB', snr1);
    text(0.5,-95,sprintf('NBW=%.0e ',NBW),'Hor','right');
    figureMagic([-0.5 0.5],1/16,4, [-100 0],10,2,[6 3],'Spectra');
    xlabel('frequency')
    ylabel('PSD (dBFS/NBW)');
    legend(msg,msg0,msg1);
    % printmif(fullfile('MIF','QMS'), [4 2], 'Helvetica8')
end

