function printmif(file,size,font,fig)
%printmif(file,size,font,fig) Print graph to an Adobe Illustrator file 
% and then use ai2mif to convert it to FrameMaker MIF format.
% ai2mif is a slightly modified version of the function of the same name
% provided by Deron Jackson <djackson@mit.edu>.

path = which(mfilename);
ai2mifpath = [ path(1:max(find(path==filesep))) 'ai2mif '];
% Handle the input arguments
parameters = {'file' 'size' 'font' 'fig'};
defaults = {'matlab' [] [] gcf };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

if ~isempty(font)
    fontsize = str2num(font(font<59));	% Numerical portion
    fontname = font(font>58);		% Alphabetical portion
    children = get(fig,'Children');
    for i = 1:length(children)
    	if isfield(children(i),'FontName');
	    set(children(i), 'FontName',fontname, 'FontSize',fontsize);
	end
    end
end

if ~isempty(size)
    set(fig,'PaperUnits','inches','PaperPosition', [0.5 0.5 0.5+size]);
end
% eval( ['print -dill -f' num2str(fig) ' ' file '.ai'] );
[path name] = fileparts(file);
%eval( ['print -dill -painters -f' num2str(fig) ' ' name '.ai'] );
print( '-dill', '-painters', ['-f' num2str(fig)], [name '.ai'] );


ai2mif(name);                           % Execute AI2MIF
delete([name '.ai']);                   % Delete the AI file
if ~isempty(path)
    MIFfile = [fullfile(cd,file) '.mif'];
    if exist(MIFfile,'file')
	delete(MIFfile);                   % Delete the old MIF file if it exists
    end
    copyfile([name '.mif'],[fullfile(cd,file) '.mif']);
    delete([name '.mif']);                 % Delete the temporary MIF file
end

