function Work = define_directories(basedir, Box_name)

Work.hmdir = basedir; % Base directory
Work.Box_base_name = Box_name; % Project name

Work.fundir = [Work.hmdir 'vbr\'];
Work.cwd=[Work.fundir '2_PLATES\4pt0_1d_plates\'];
Work.veldir =  [Work.fundir '6_FitVobs\Functions_Emily\vel_models\'];
Work.labdir =  [Work.fundir '6_FitVobs\Functions_Emily\LAB_models\'];
Work.valsdir = [Work.fundir '6_FitVobs\Functions_Emily\hardwired_values\'];
addpath(Work.cwd);
addpath([Work.fundir '6_FitVobs/pyFits_v0p1/']);

end