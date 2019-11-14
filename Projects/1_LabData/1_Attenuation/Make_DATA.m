clear Data
% make the data structure with all experimental results !
% The Cooper one is Gribb and Cooper, 1998 with two temps (three points each).
% The Jackson et al., 2002 data has 10 points and the Tan et al., 2001 has 7.
% Sundberg data is two temps but I remove the two or three high freq. points
% that make the peak (this is probably not a legit thing to do, which is why
% I don't use Sundberg data in the JGR paper.
% But his data make a nice low freq. tail to compare to the andrade curve
% so you should use it). columns that say either f/fm or normalized f are
% the ones I use for the master curve. in all cases I used the viscosities
% reported in those papers.

% November 8, 2019 (Ben):
% the data is stored outside of the vbr repository because we do not want to release
% the data with the codes without explicit permission... so it must be read in..
% [AND IT MUST BE MOVED OUT OF THIS CODE !! ]
% re

% ======================================================
% FAUL AND JACKSON 2015 DATA (that includes previous papers, right?)
datadir = '../../../../vbrWork/expt_data/3_attenuation/' ;
FJ_datadir = [datadir,'FJ2015_data/'] ;
addpath(FJ_datadir) ;
G_filelist = {"T900_SanCarlos_G.csv", "T1000_SanCarlos_G.csv", "T1100_SanCarlos_G.csv", "T1200_SanCarlos_G.csv"}
Qinv_filelist = {"T900_SanCarlos_Qinv.csv","T1000_SanCarlos_Qinv.csv","T1100_SanCarlos_Qinv.csv","T1200_SanCarlos_Qinv.csv"};
T_C_vec = [900 1000 1100 1200]; % ideally read these in from the filenames :)


for iT=1: length(T_C_vec)

  T_C = T_C_vec(iT)
  Gfilename = G_filelist{iT} ;
  disp(Gfilename);
  data_G = load(Gfilename);
  logPer = transpose(data_G(:,1)) ;
  f = 1./(10.^logPer);
  Qfilename = Qinv_filelist{iT};
  disp(Qfilename);
  data_Qinv = load(Qfilename);

  Data.FaulJax15(iT).exptCond.T_C = T_C ;
  Data.FaulJax15(iT).exptCond.P_GPa = 0.300 ; % confining pressure
  Data.FaulJax15(iT).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
  Data.FaulJax15(iT).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
  Data.FaulJax15(iT).exptCond.dg_0 = 17.1 ; % grain size, in microns
  Data.FaulJax15(iT).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
  Data.FaulJax15(iT).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
  Data.FaulJax15(iT).exptCond.phi_0 = 0.00 ; % melt fraction
  Data.FaulJax15(iT).exptCond.logPer = logPer ; %
  Data.FaulJax15(iT).exptCond.f = f ;
  Data.FaulJax15(iT).exptCond.logf = log10(f) ;
  Data.FaulJax15(iT).Results.log10_Qinv = data_Qinv(:,2) ;
  Data.FaulJax15(iT).Results.G = data_G(:,2) ;

end



% ======================================================
% BUNTON and COOPER DATA

datadir = '../../../../vbrWork/expt_data/3_attenuation/' ;
BC_datadir = [datadir,'Bunton/'] ;
addpath(BC_datadir) ;
%
% Qinv_freq_1200C.csv  Qinv_freq_1250C.csv   Qinv_freq_1300C.csv
%G_filelist = {'logflogQinv_d4_1200C.csv', 'G_freq_1250C.csv', 'G_freq_1300C.csv'} ;
Qinv_filelist = {'logflogQinv_d4_1200C.csv', 'logflogQinv_d4_1225C.csv', 'logflogQinv_d4_1250C.csv', 'logflogQinv_d4_1285C.csv', 'logflogQinv_d4_1300C.csv'};
T_C_vec = [1200 1225 1250 1285 1300]; % ideally read these in from the filenames :)


for iT=1:length(T_C_vec)

  T_C = T_C_vec(iT) ;
  %Gfilename = G_filelist{iT} ;
  %disp(Gfilename);
  %data_G = load(Gfilename);
  %f = transpose(data_G(:,1)) ;
  %Per = 1./f;
  %logPer = log10(Per) ;

  Qfilename = Qinv_filelist{iT};
  disp(Qfilename);
  data_Qinv = load(Qfilename);
  log10f = transpose(data_Qinv(:,1)) ;
  f = 10.^log10f ;
  Per = 1./f;
  logPer = log10(Per) ;

  Data.BuntonCoop01(iT).Results.log10_Qinv = data_Qinv(:,2) ;
  %Data.BuntonCoop10(iT).Results.G = data_G(:,2) ;

  Data.BuntonCoop01(iT).exptCond.T_C = T_C ;
  Data.BuntonCoop01(iT).exptCond.P_GPa = 101325*1e-9 ; % confining pressure-- Room pressure
  Data.BuntonCoop01(iT).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
  Data.BuntonCoop01(iT).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
  Data.BuntonCoop01(iT).exptCond.dg_0 = 4.0 ; % grain size, in microns
  Data.BuntonCoop01(iT).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
  Data.BuntonCoop01(iT).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
  Data.BuntonCoop01(iT).exptCond.phi_0 = 0.00 ; % melt fraction
  Data.BuntonCoop01(iT).exptCond.logPer = logPer ; %
  Data.BuntonCoop01(iT).exptCond.f = f ;
  Data.BuntonCoop01(iT).exptCond.logf = log10f ;

end


% ======================================================
% YAMAUCHI-TAKEI 2016 data

datadir = '../../../../vbrWork/expt_data/3_attenuation/' ;
YT_datadir = [datadir,'Yamauchi2016/'] ;
addpath(YT_datadir) ;

%
% Qinv_freq_1200C.csv  Qinv_freq_1250C.csv   Qinv_freq_1300C.csv
%G_filelist = {'logflogQinv_d4_1200C.csv', 'G_freq_1250C.csv', 'G_freq_1300C.csv'} ;
%Qinv_filelist = {'logflogQinv_d4_1200C.csv', 'logflogQinv_d4_1225C.csv', 'logflogQinv_d4_1250C.csv', 'logflogQinv_d4_1285C.csv', 'logflogQinv_d4_1300C.csv'};
Qinv_filelist = dir([YT_datadir,'YT16_s40_fQinv_*']) ;
G_filelist = dir([YT_datadir,'YT16_s40_fE_*']) ;

T_C_vec = [8 13 18 23 29 35 39 47] ; % ideally read these in from the filenames :)


for iT=1:length(T_C_vec)

  T_C = T_C_vec(iT) ;

  Gfilename = G_filelist(iT).name ;
  disp(Gfilename);
  data_G = load(Gfilename);
  f = transpose(data_G(:,1)) ;
  Per = 1./f;
  logPer = log10(Per) ;

  Qfilename = Qinv_filelist(iT).name;
  disp(Qfilename);
  data_Qinv = load(Qfilename);
  f = transpose(data_Qinv(:,1)) ;
  Per = 1./f;
  logPer = log10(Per) ;

  Data.YT16(iT).Results.Qinv = data_Qinv(:,2) ;
  Data.YT16(iT).Results.log10_Qinv = log10(Data.YT16(iT).Results.Qinv) ;
  Data.YT16(iT).Results.G = data_G(:,2) ;

  Data.YT16(iT).exptCond.T_C = T_C ;
  Data.YT16(iT).exptCond.P_GPa = 101325*1e-9 ; % confining pressure-- Room pressure
  Data.YT16(iT).exptCond.Gu = 2.6 ; % reference G ! not sure what this is.. placeholder
  Data.YT16(iT).exptCond.rho = 1.011e3 ; % reference density ! not sure what this is.. placeholder
  Data.YT16(iT).exptCond.dg_0 = 24.4 ; % grain size, in microns
  Data.YT16(iT).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
  Data.YT16(iT).exptCond.sig_0 = 1e3 ; % stress,  Pa (cnvrted to MPa in spine)
  Data.YT16(iT).exptCond.phi_0 = 0.00 ; % melt fraction
  Data.YT16(iT).exptCond.logPer = logPer ; %
  Data.YT16(iT).exptCond.f = f ;
  Data.YT16(iT).exptCond.logf = log10f ;

end

% ======================================================
% McCARTHY_etal-TAKEI 2011 data

datadir = '../../../../vbrWork/expt_data/3_attenuation/' ;
McCT_datadir = [datadir,'McCT11/'] ;
addpath(McCT_datadir) ;

Qinv_filelist = dir([McCT_datadir,'McCT11_fQinv_*']) ;
%G_filelist = dir([McCT_datadir,'YT16_s40_fE_*']) ;

T_C_vec = [23.0] ; % ideally read these in from the filenames :)
d_vec = [4.3 6.3] ;

for iT=1:length(d_vec)

  T_C = T_C_vec(1) ;
  d = d_vec(iT) ;

  % Gfilename = G_filelist(iT).name ;
  % disp(Gfilename);
  % data_G = load(Gfilename);
  % f = transpose(data_G(:,1)) ;
  % Per = 1./f;
  % logPer = log10(Per) ;

  Qfilename = Qinv_filelist(iT).name;
  disp(Qfilename);
  data_Qinv = load(Qfilename);
  f = transpose(data_Qinv(:,1)) ;
  Per = 1./f;
  logPer = log10(Per) ;

  Data.McCT11(iT).Results.Qinv = data_Qinv(:,2) ;
  Data.McCT11(iT).Results.log10_Qinv = log10(Data.McCT11(iT).Results.Qinv) ;
  %Data.YT16(iT).Results.G = data_G(:,2) ;

  Data.McCT11(iT).exptCond.T_C = T_C ;
  Data.McCT11(iT).exptCond.P_GPa = 101325*1e-9 ; % confining pressure-- Room pressure
  Data.McCT11(iT).exptCond.Gu = 2.6 ; % reference G ! not sure what this is.. placeholder
  Data.McCT11(iT).exptCond.rho = 1.011e3 ; % reference density ! not sure what this is.. placeholder
  Data.McCT11(iT).exptCond.dg = d_vec(iT) ; % grain size, in microns
  Data.McCT11(iT).exptCond.dg_0 = d_vec(1) ; % grain size, in microns
  Data.McCT11(iT).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
  Data.McCT11(iT).exptCond.sig_0 = 1e3 ; % stress,  Pa (cnvrted to MPa in spine)
  Data.McCT11(iT).exptCond.phi_0 = 0.00 ; % melt fraction
  Data.McCT11(iT).exptCond.logPer = logPer ; %
  Data.McCT11(iT).exptCond.f = f ;
  Data.McCT11(iT).exptCond.logf = log10f ;

end




% ===================
save('./ExptData.mat', 'Data')
% ===================


% ======================================================
return
% ======================================================


% ======================================================
% SUNDBERG and COOPER DATA

datadir = '../../../../vbrWork/expt_data/3_attenuation/' ;
SC_datadir = [datadir,'SundCoop2010_data/'] ;
addpath(SC_datadir) ;
%
% Qinv_freq_1200C.csv  Qinv_freq_1250C.csv   Qinv_freq_1300C.csv
G_filelist = {'G_freq_1200C.csv', 'G_freq_1250C.csv', 'G_freq_1300C.csv'} ;
Qinv_filelist = {'Qinv_freq_1200C.csv', 'Qinv_freq_1250C.csv', 'Qinv_freq_1300C.csv'};
T_C_vec = [1200 1250 1300]; % ideally read these in from the filenames :)


for iT=1:length(T_C_vec)

  T_C = T_C_vec(iT)
  Gfilename = G_filelist{iT} ;
  disp(Gfilename);
  data_G = load(Gfilename);
  f = transpose(data_G(:,1)) ;
  Per = 1./f;
  logPer = log10(Per) ;
  Qfilename = Qinv_filelist{iT};
  disp(Qfilename);
  data_Qinv = load(Qfilename);

  Data.SundCoop10(iT).Results.log10_Qinv = data_Qinv(:,2) ;
  Data.SundCoop10(iT).Results.G = data_G(:,2) ;

  Data.SundCoop10(iT).exptCond.T_C = T_C ;
  Data.SundCoop10(iT).exptCond.P_GPa = 0.300 ; % confining pressure
  Data.SundCoop10(iT).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
  Data.SundCoop10(iT).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
  Data.SundCoop10(iT).exptCond.dg_0 = 17.1 ; % grain size, in microns
  Data.SundCoop10(iT).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
  Data.SundCoop10(iT).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
  Data.SundCoop10(iT).exptCond.phi_0 = 0.00 ; % melt fraction
  Data.SundCoop10(iT).exptCond.logPer = logPer ; %
  Data.SundCoop10(iT).exptCond.f = f ;
  Data.SundCoop10(iT).exptCond.logf = log10(f) ;

end

save('./ExptData.mat', 'Data')




% ========
return
% ========
% ALL THE REST MUST BE REPLACED AND DELETED !



%% TAN & FAUL & JACKSON ()
% add in the rest of this data !

Data.TanJax(1).exptCond.T_C = 1300 ;
Data.TanJax(1).exptCond.P_GPa = 0.300 ; % confining pressure
Data.TanJax(1).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.TanJax(1).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.TanJax(1).exptCond.dg_0 = 10. ; % CHECK THIS ! ; % grain size, in microns
Data.TanJax(1).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.TanJax(1).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.TanJax(1).exptCond.phi_0 = 0.00 ; % melt fraction
Data.TanJax(1).exptCond.f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ];
Data.TanJax(1).exptCond.logf = log10(Data.TanJax(1).exptCond.f) ; % log10(Data.SundCoop(2).exptCond.f) ;
Data.TanJax(1).Results.Qinv = [ 1.7473 1.2105 0.9289 0.7106 0.5611 0.4459 0.3586 0.3251 0.2992 0.3001 ] ;
Data.TanJax(1).Results.G = [ 4.9 7.4 9.8 12.5 15.0 17.4 20.0 22.3 23.8 26.0 ] ;




%% GRIBB AND COOPER
% T_params = [ 1225 1285 ] + 273 ;
% BH, June 7, 2016: not sure what this is ?
% see below, double check the temperatures.

Data.GribCoop(1).exptCond.T_C = 1200 ; % or 1225 as above?
Data.GribCoop(1).exptCond.P_GPa = 0.001 ; % confining pressure
Data.GribCoop(1).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.GribCoop(1).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.GribCoop(1).exptCond.dg_0 = 10 ; % actual=2.8 ; % grain size, in microns
Data.GribCoop(1).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.GribCoop(1).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.GribCoop(1).exptCond.phi_0 = 0.00 ; % melt fraction
Data.GribCoop(1).exptCond.logf = [-1.8789 -1.2063 -0.89238 ] ;
Data.GribCoop(1).exptCond.f = [ 0.013215 0.062190 0.12812 ] ;
Data.GribCoop(1).Results.G = [ 24.211 29.474 31.579 ] ;
Data.GribCoop(1).Results.logQinv = [ -0.50000 -0.72000 -0.83333 ] ;
Data.GribCoop(1).Results.Qinv = [ 0.31623 0.19055 0.14678 ] ;

Data.GribCoop(2).exptCond.T_C = 1300 ; % or 1285 as above
Data.GribCoop(2).exptCond.P_GPa = 0.001 ; % confining pressure
Data.GribCoop(2).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.GribCoop(2).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.GribCoop(2).exptCond.dg_0 = 10 ; % actual=2.8 ; % grain size, in microns
Data.GribCoop(2).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.GribCoop(2).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.GribCoop(2).exptCond.phi_0 = 0.00 ; % melt fraction
Data.GribCoop(2).exptCond.logf = [ -2.2377 -1.8789 -0.4888 ] ;
Data.GribCoop(2).exptCond.f = [ 0.005785 0.013215 0.32450 ] ;
Data.GribCoop(2).Results.G = [ 11.579 15.789 26.842 ] ;
Data.GribCoop(2).Results.logQinv = [ 0.0000 -0.20000 -0.76667 ] ;
Data.GribCoop(2).Results.Qinv = [ 1.0000 0.63096 0.17113 ] ;

%% SUNDBERG AND COOPER, 2010

Data.SundCoop(1).exptCond.T_C = 1200 ;
Data.SundCoop(1).exptCond.P_GPa = 0.001 ; % confining pressure
Data.SundCoop(1).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.SundCoop(1).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.SundCoop(1).exptCond.dg_0 = 2.8 ; % actual=2.8 ; % grain size, in microns
Data.SundCoop(1).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.SundCoop(1).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.SundCoop(1).exptCond.phi_0 = 0.00 ; % melt fraction
Data.SundCoop(1).exptCond.f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ];
Data.SundCoop(1).exptCond.logf = log10(Data.SundCoop(1).exptCond.f) ;
Data.SundCoop(1).Results.Qinv = [ 0.4301 0.35230 0.28160 0.21830 0.19740 0.16490 0.14690 0.17370 0.18610 0.22820 ] ;
Data.SundCoop(1).Results.G = [21.2 23.8 25.3 27.3 29.0 30.6 32.0 32.9 34.6 35.9 ] ;

Data.SundCoop(2).exptCond.T_C = 1250 ;
Data.SundCoop(2).exptCond.P_GPa = 0.001 ; % confining pressure
Data.SundCoop(2).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.SundCoop(2).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.SundCoop(2).exptCond.dg_0 = 2.8 ; % actual=2.8 ; % grain size, in microns
Data.SundCoop(2).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.SundCoop(2).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.SundCoop(2).exptCond.phi_0 = 0.00 ; % melt fraction
Data.SundCoop(2).exptCond.f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ];
Data.SundCoop(2).exptCond.logf = log10(Data.SundCoop(2).exptCond.f) ;
Data.SundCoop(2).Results.Qinv = [ 0.816 0.6354 0.4964 0.4101 0.3374 0.2791 0.2333 0.2124 0.2175 0.2405] ;
Data.SundCoop(2).Results.G = [ 12.2 14.9 17.8 20.0 22.4 24.7 26.6 28.3 29.9 31.1 ] ;

Data.SundCoop(3).exptCond.T_C = 1300 ;
Data.SundCoop(3).exptCond.P_GPa = 0.001 ; % confining pressure
Data.SundCoop(3).exptCond.Gu = 65 ; % reference G ! not sure what this is.. placeholder
Data.SundCoop(3).exptCond.rho = 3300 ; % reference density ! not sure what this is.. placeholder
Data.SundCoop(3).exptCond.dg_0 = 2.8 ; % actual=2.8 ; % grain size, in microns
Data.SundCoop(3).exptCond.Ch2o_0 = 0 ; % water content (wt %?)
Data.SundCoop(3).exptCond.sig_0 = 1e5 ; % stress,  Pa (cnvrted to MPa in spine)
Data.SundCoop(3).exptCond.phi_0 = 0.00 ; % melt fraction
Data.SundCoop(3).exptCond.f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ];
Data.SundCoop(3).exptCond.logf = log10(Data.SundCoop(2).exptCond.f) ;
Data.SundCoop(3).Results.Qinv = [ 1.7473 1.2105 0.9289 0.7106 0.5611 0.4459 0.3586 0.3251 0.2992 0.3001 ] ;
Data.SundCoop(3).Results.G = [ 4.9 7.4 9.8 12.5 15.0 17.4 20.0 22.3 23.8 26.0 ] ;
