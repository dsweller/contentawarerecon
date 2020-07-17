function espiritpath = setup_ESPIRiT(espiritpath)
% setup ESPIRiT toolbox (download from web -- http://people.eecs.berkeley.edu/~mlustig/Software.html)

zpad_path = which('zpad','-all');
[pathstrs,filenames,fileexts] = cellfun(@(p) fileparts(p), zpad_path,'UniformOutput',false);
goodfiles = strcmp(filenames,'zpad') & strcmp(fileexts,'.m') & ~cellfun(@isempty,pathstrs);
pathstrs = pathstrs(goodfiles);
pathparts = cellfun(@(p) strsplit(p,filesep),pathstrs,'UniformOutput',false);
goodpaths = cellfun(@(p) strcmp(p{end},'utils'),pathparts);
if any(goodpaths)
    igoodpath = find(goodpaths,1,'first');
    pathparts = pathparts{igoodpath};
    espiritpath = fullfile(pathparts{1:end-1});
else
    % need to load ESPIRiT toolbox files (not full toolbox)
    if ~exist('espiritpath','var') || isempty(espiritpath), espiritpath = ['..' filesep 'ESPIRiT']; end
    while exist(espiritpath,'dir') ~= 7 || exist([espiritpath filesep 'setPath.m'],'file') ~= 2
        if usejava('awt')
            espiritpath = uigetdir('.','Locate ESPIRiT Toolbox...');
        else
            espiritpath = input('Path to ESPIRiT: ','s');
        end
        if isequal(espiritpath,0)
            return; % later function calls may fail!
        end
    end
    oldp = cd(espiritpath);
	addpath(espiritpath);
    addpath(fullfile(espiritpath,'utils'));
    addpath(fullfile(espiritpath,'SPIRiT_code'));
    addpath(fullfile(espiritpath,'SAKE_code'));
    addpath(fullfile(espiritpath,'ESPIRiT_code'));
    addpath(fullfile(espiritpath,'nufft_files'));
    cd(oldp);
end

end
