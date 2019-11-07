%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mantle_extrap.m
%
% - generates a range of mantle conditions from thermal modesl
% - calculates mechnaical properties from the thermal models
% - compares mechanical properties between methods for mantle conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% put VBR in the path
path_to_top_level_vbr='../../';
addpath(path_to_top_level_vbr)
vbr_init

% add the local dirs
addpath('./data')
addpath('./functions')

% generate mantle conditions from thermal model (or load if it exists)
ThermalSettings.Tpots=1200:100:1500;
Box = genPullThermalModels(fullfile(pwd,'data/thermalModelBox.mat'),ThermalSettings);

% run Box through VBR calculator (or load if it exists)
VBRsettings.ane_meths={'AndradePsP';'MTH2011';'eBurgers';'YT2016_solidus'};
VBRsettings.freqs=[0.01, 0.1];
VBRsettings.phi0=0.01; % phi when T > Tsol
Box = genPullVBRdata(Box,fullfile(pwd,'data/VBR_Box.mat'),VBRsettings);
