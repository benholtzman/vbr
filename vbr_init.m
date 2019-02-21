function vbr_init(VBR_version,PLATES_version)
% adds all relevant VBR paths to the matlab path
%
% if within the vbr top level directory, just call:
%   vbr_init
%
% if elsewhere (or within scripts), you need to add the path to the top level
% vbr directory first:
%   addpath('/path/to/vbr')
%   vbr_init 
%
% optional inputs:
%     VBR_version  string, if not in VBR_versions, will use default
%     PLATES_version  string, if not in PLATES_versions, will use default

% define available versions, defaults
  VBR_versions={'VBR_v0p95','VBR_stable'};
  VBR_default=VBR_versions{1};
  PLATES_versions={'4pt0_1d_plates'};
  PLATE_default=PLATES_versions{1};

% get full path to vbr, regardless of where vbr_init is called from
  p=mfilename('fullpath'); % full path of vbr_init without extension
  [filepath,name,ext] = fileparts([p,'.m']);
  vbr_dir=filepath; % remove filename from fullpath

% put VBR core in the path
  VBR_dir=VBR_default;
  if nargin >= 1
    if any(strcmp(VBR_versions,VBR_version))
      VBR_dir=VBR_version;
    end
  end
  disp(['initializing with VBR_version: ',VBR_dir])
  VBR_path=fullfile(vbr_dir,'vbr','4_VBR',VBR_dir);
  addpath(genpath(VBR_path));

% put plates in path
  plate_dir=PLATE_default;
  if nargin > 1
    if any(strcmp(PLATES_versions,PLATES_version))
      plate_dir=PLATES_version;
    end
  end
  plate_path=fullfile(vbr_dir,'vbr','2_PLATES',plate_dir);
  addpath(genpath(plate_path));

% put remainder of folders in path
  fo2add={'0_COOKBOOK','1_LabData','6_FitVobs'};
  for i_fo = 1:numel(fo2add)
      fo=fo2add{i_fo};
      path2add=fullfile(vbr_dir,'vbr',fo);
      addpath(genpath(path2add));
  end

  disp('VBR calculator initialized');
end
