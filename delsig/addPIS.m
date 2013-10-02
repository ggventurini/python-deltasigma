function addPIS
% Add the PosInvSet subdirectory of the toolbox to the current path
tb_dir = which('dsdemo5');
slash = find(tb_dir==filesep);
if isempty(slash)
    tb_dir = '';
else
    tb_dir = tb_dir(1:slash(end));
end
addpath([tb_dir 'PosInvSet']);
