function vbr_init(varargin)
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
% optional keyword-value pair inputs:
%   'VBR_version','VBR_v0p95'  to use a VBR version other than the default
%   'PLATES_version','4pt0_1d_plates'  to use a thermal model other than the default
%

% define available versions, defaults
  ValidOpts=struct();
  ValidOpts.vbr_version={'VBR_v0p95','VBR_stable'};
  ValidOpts.plates_version={'4pt0_1d_plates'};
  Options=struct('vbr_version',ValidOpts.vbr_version{1},...
                 'plates_version',ValidOpts.plates_version{1});

% get full path to vbr, regardless of where vbr_init is called from
  p=mfilename('fullpath'); % full path of vbr_init without extension
  [filepath,name,ext] = fileparts([p,'.m']);
  vbr_dir=filepath; % remove filename from fullpath

% add the vbr/support directory and validate input options
  addpath(genpath(fullfile(vbr_dir,'vbr','support')));
  Options=validateStructOpts('vbr_init',varargin,Options,ValidOpts);
  disp(['initializing with vbr_version: ',Options.vbr_version])

% collect all the subdirectories under ./vbr/ to add
  subDirs2add={fullfile('4_VBR',Options.vbr_version),...
               fullfile('2_PLATES',Options.plates_version),...
               '0_COOKBOOK','1_LabData',...
               fullfile('6_FitVobs','Functions_Bayesian'),...
               fullfile('6_FitVobs','pyFits_v0p1')};
  for i_fo = 1:numel(subDirs2add)
     fo=subDirs2add{i_fo};
     path2add=fullfile(vbr_dir,'vbr',fo);
     if exist(path2add,'dir')
       addpath(genpath(path2add));
     end
  end

  disp('VBR calculator initialized');
end
