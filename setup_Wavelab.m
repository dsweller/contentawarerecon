function wavelabpath = setup_Wavelab(wavelabpath)
% setup Wavelab toolbox v850 (download from web -- http://www-stat.stanford.edu/~wavelab)

wavepath_path = which('WavePath','-all');
[pathstrs,filenames,fileexts] = cellfun(@(p) fileparts(p), wavepath_path,'UniformOutput',false);
goodfiles = strcmp(filenames,'WavePath') & strcmp(fileexts,'.m') & ~cellfun(@isempty,pathstrs);
pathstrs = pathstrs(goodfiles);
if any(goodfiles)
    igoodfile = find(goodfiles,1,'first');
    wavelabpath = pathstrs{igoodfile};
else
    % need to load WaveLab toolbox files (not full toolbox)
    if ~exist('wavelabpath','var') || isempty(wavelabpath), wavelabpath = ['..' filesep 'Wavelab850']; end
    while exist(wavelabpath,'dir') ~= 7 || exist([wavelabpath filesep 'WavePath.m'],'file') ~= 2
        if usejava('awt')
            wavelabpath = uigetdir('.','Locate Wavelab Toolbox...');
        else
            wavelabpath = input('Path to Wavelab: ','s');
        end
        if isequal(wavelabpath,0)
            return; % later function calls may fail!
        end
    end
    oldp = cd(wavelabpath);
	addpath(wavelabpath);
    addpath(fullfile(wavelabpath,'Orthogonal'));
    cd(oldp);
end

end
