function changeFig(fontsize,linewidth,markersize)
% changeFig(fontsize,linewidth,markersize)	Change the settings
% for all objects in the current figure window
% Handle the input arguments
parameters = {'fontsize' 'linewidth' 'markersize'};
defaults = { 9 1 6 };
for i=1:length(defaults)
    parameter = char(parameters(i));
    if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
     eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
        eval([parameter '=defaults{i};'])
    end
end

children = get(gcf,'children');
children = children(:)';
for child=children
    type = get(child,'type');
    if type=='axes'
	set(child,'fontsize',fontsize)
	axisChildren = get(child,{'children' 'xlabel' 'ylabel' 'title'});
	axisChildren = [axisChildren{1}' axisChildren{2:4}];
	for axisChild=axisChildren
	    type = get(axisChild,'type');
	    switch type
		case 'line'
		    set(axisChild,'linewidth',linewidth);
		    set(axisChild,'markersize',markersize);
		case 'text'
		    set(axisChild,'fontsize',fontsize);
		otherwise
		    ;
	    end
	end
    else
	fprintf(1,'type is %s (not ''axes'')\n',type);
    end
end
