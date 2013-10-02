function axis_handle = DocumentNTF(arg1,osr,f0,quadrature)
%axis_handle = DocumentNTF(ntf|ABCD|mod_struct,osr=64,f0=0,quadrature=0)   Plot the NTF's poles and zeros as well as its frequency-response
%
% The first argument is either the NTF, ABCD matrix or a struct containing 
% ntf, osr=64, f0=0, quadrature=0 and optionally stf. 
% If the first argument is a struct, then no other arguments should be supplied.
% If the first argument is ABCD, the stf is also plotted.

% Handle the input arguments
parameters = {'arg1' 'osr' 'f0' 'quadrature'};
defaults = { [] 64 0 0 };
for arg_ii=1:length(defaults)
    parameter = parameters{arg_ii};
    if arg_ii>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{arg_ii};'])
    end
end
if nargin==1 & isstruct(arg1)
    flattenStruct(arg1);
    for arg_ii=2:length(defaults)
        parameter = parameters{arg_ii};
        if ~exist(parameter,'var')
            eval([parameter '=defaults{arg_ii};'])
        end
    end
elseif isnumeric(arg1)
    ABCD = arg1;
    [ntf stf] = calculateTF(ABCD);
else
    ntf = arg1;
end
set(gcf,'name','NTF')
set(gcf,'numbertitle','off');
set(gcf,'MenuBar','none');

% Plot poles and zeros
clf;
ax_handle(1) = subplot('position',[.05 .15 .25 .65]);
plotPZ(ntf,'b',6);
axis([-1.1 1.1 -1.1 1.1]);
set(gca,'XTick', [-1:0.5:1]);
set(gca,'XTickLabel', {'-1','','0','','1'});
set(gca,'YTick', [-1:0.5:1]);
set(gca,'YTickLabel', {'-1','','0','','1'});
title('Poles and Zeros');
% Frequency response
ax_handle(2) = subplot('position', [0.37,0.1,0.6 0.8]);
f = ds_freq(osr,f0,quadrature);
z = exp(2i*pi*f);
H = dbv(evalTF(ntf,z));
plot(f,H,'b');
if exist('stf','var') & ~isempty(stf)
    set(gcf,'name','NTF and STF')
    G = dbv(evalTF(stf,z));
    hold on;
    plot(f,G,'m');
    hold off;
end
[f1 f2] = ds_f1f2(osr,f0,quadrature);
NG0 = dbv(rmsGain(ntf,f1,f2));
hold on;
plot([f1 f2],NG0*[1 1],'k','Linewidth', 3);
if f0==0
    text(0.5/osr,NG0,sprintf('  %.0fdB',NG0),'Vert','Mid','Hor','Left');
else
    text(f0,NG0+1,sprintf('%.0fdB',NG0),'Vert','Bot','Hor','Cen');
end
msg = sprintf(' Inf-norm of H=%.2f\n 2-norm of H=%.2f', infnorm(ntf), rmsGain(ntf,0,1));
if f0<0.25
    text(0.48,0,msg,'Vert','Top','Hor','Right');
else
    text(f_left,0,msg,'Vert','Top');
end
if quadrature
    ING0 = dbv(rmsGain(ntf,-f1,-f2));
    plot(-[f1 f2],ING0*[1 1],'k','Linewidth', 3);
    text(-f0,ING0+1,sprintf('%.0fdB',ING0),'Vert','Bot','Hor','Cen');
    f_left = -0.5;
else
    f_left = 0;
end
figureMagic([f_left,0.5],1/16,2, [-80 15],10,2);
xlabel('frequency');
% text(0.25,-85,'Normalized frequency (1\rightarrow f_s)','hor','center');
title('Frequency Response');

if nargout>0
    axis_handle = ax_handle;
end
