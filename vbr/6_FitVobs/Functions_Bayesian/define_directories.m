function Work = define_directories(basedir, Box_name)

Work.hmdir = basedir; % Base directory
Work.Box_base_name = Box_name; % Project name

Work.fundir = [Work.hmdir 'vbr/vbr/'];
Work.cwd=[Work.fundir '2_PLATES/4pt0_1d_plates/'];
Work.veldir =  [Work.hmdir 'Seismic_Inputs/vel_models/'];
Work.labdir =  [Work.hmdir 'Seismic_Inputs/LAB_models/'];
Work.valsdir = [Work.hmdir 'Seismic_Inputs/hardwired_values/'];
addpath(Work.cwd);
addpath([Work.fundir '6_FitVobs/pyFits_v0p1/']);

end