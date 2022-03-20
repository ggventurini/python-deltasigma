function axis_handle = DocumentNTF(arg1,osr,f0,quadrature,sizes,plotFreqRespOnly)
%axis_handle = DocumentNTF(ntf|ABCD|mod_struct,osr=64,f0=0,quadrature=0,sizes,plotFreqRespOnly=1)   
%
% The first argument is either the NTF, ABCD matrix or a struct containing 
% ntf, osr=64, f0=0, quadrature=0 and optionally stf. 
% If the first argument is a struct, then no other arguments should be supplied.
% If the first argument is ABCD, the stf is also plotted.

% Handle the input arguments
parameters = {'arg1' 'osr' 'f0' 'quadrature' 'sizes' 'plotFreqRespOnly'};
defaults = { [] 64 0 0 NaN 0};
for arg_ii=1:length(defaults)
    parameter = parameters{arg_ii};
    if arg_ii>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) ) %#ok<OR2,AND2>
        eval([parameter '=defaults{arg_ii};'])
    end
end
if nargin==1 && isstruct(arg1)
    for arg_ii=2:length(defaults)
        parameter = parameters{arg_ii};
        if isfield(arg1,parameter)
            eval([parameter '= arg1.(parameter);'])
	else
            eval([parameter '=defaults{arg_ii};'])
        end
    end
elseif isnumeric(arg1)
    ABCD = arg1;
    [ntf, stf] = calculateTF(ABCD);
else
    ntf = arg1;
end
if isnumeric(sizes) && all(isnan(sizes))
    clear sizes;
    sizes.lw = 1;        % LineWidth
    sizes.ms = 5;        % MarkerSize
    sizes.fs = 12;       % FontSize
    sizes.fw = 'normal'; % FontWeight
end
logplot = f0==0;
if quadrature
    f_left = -0.5;
elseif logplot
    f_left = 1e-3;
else
    f_left = 0;
end
    
set(gcf,'name','NTF')
set(gcf,'numbertitle','off');
set(gcf,'MenuBar','none');
clf;
if ~plotFreqRespOnly
    % Plot poles and zeros
    clf;
    ax_handle(1) = subplot('position',[.05 .15 .25 .65]);
    plotPZ(ntf,'b',6);
    axis([-1.1 1.1 -1.1 1.1]);
    set(gca,'XTick', (-1:0.5:1));
    set(gca,'XTickLabel', {'-1','','0','','1'});
    set(gca,'YTick', (-1:0.5:1));
    set(gca,'YTickLabel', {'-1','','0','','1'});
    title('Poles and Zeros');
    ax_handle(2) = subplot('position', [0.37 0.12 0.6 0.78]);
else
    ax_handle(2) = subplot('position', [0.1 0.15 0.85 0.75]);
end
% Frequency response
f = ds_freq(osr,f0,quadrature);
z = exp(2i*pi*f);
H = dbv(evalTF(ntf,z));
if logplot
    semilogx(f,H,'b','Linewidth',sizes.lw);
else
    plot(f,H,'b','Linewidth',sizes.lw);
end
if exist('stf','var') && ~isempty(stf)
    set(gcf,'name','NTF and STF')
    G = dbv(evalTF(stf,z));
    hold on;
    plot(f,G,'m','Linewidth',sizes.lw);
    hold off;
end
[f1,f2] = ds_f1f2(osr,f0,quadrature);
NG0 = dbv(rmsGain(ntf,f1,f2));
hold on;
plot([max(f1,f_left) f2],NG0*[1 1],'--k','Linewidth', sizes.lw+1 );
if f0==0 & logplot
    h = text(sqrt(f_left*0.5/osr),NG0+2,sprintf('  %.0fdB',NG0),'Vert','Bot','Hor','Cen');
else
    h = text(f2,NG0-1,sprintf(' %.0fdB',NG0),'Vert','Mid','Hor','Left');
end
set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
msg = sprintf(' Inf-norm of H = %.2f\n 2-norm of H = %.2f', infnorm(ntf), rmsGain(ntf,0,1));
if f0<0.25
    h = text(0.48,0,msg,'Vert','Top','Hor','Right');
else
    h = text(f_left,0,msg,'Vert','Top');
end
set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
if quadrature
    ING0 = dbv(rmsGain(ntf,-f1,-f2));
    plot(-[f1 f2],ING0*[1 1],'k','Linewidth', 3);
    h = text(-f0,ING0+1,sprintf('%.0fdB',ING0),'Vert','Bot','Hor','Cen');
end
set(h, 'FontSize',sizes.fs, 'FontWeight',sizes.fw);
if logplot
    figureMagic([f_left,0.5],[],[], [-100 15],10,2);
else
    figureMagic([f_left,0.5],1/16,2, [-100 15],10,2);
end
xlabel('Normalized Frequency');
% text(0.25,-85,'Normalized frequency (1\rightarrow f_s)','hor','center');
title('Frequency Response');

if nargout>0
    axis_handle = ax_handle;
end
