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
addpath('./datasets/FJ2015_data/')
G_filelist = {"T900_SanCarlos_G.csv", "T1000_SanCarlos_G.csv", "T1100_SanCarlos_G.csv", "T1200_SanCarlos_G.csv"}
Qinv_filelist = {"T900_SanCarlos_Qinv.csv","T1000_SanCarlos_Qinv.csv","T1100_SanCarlos_Qinv.csv","T1200_SanCarlos_Qinv.csv"};
T_C_vec = [900 1000 1100 1200];


for iT=1:4

  T_C = T_C_vec(iT)
  Gfilename = string(G_filelist(iT)) ; % in Octave, use 'char' instead of 'string'.
  disp(Gfilename);
  data_G = load(Gfilename);

  Qfilename = string(Qinv_filelist(iT)); % in Octave, use 'char' instead of 'string'.
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
  %Data.FaulJax15(iT).exptCond.f = [ 0.0056000 0.0100000 0.017800 0.031600 0.056200 0.10000 0.17780 0.31600 0.56230 1.0000 ];
  logPer = transpose(data_G(:,1)) ;
  f = 1./(10.^logPer);
  Data.FaulJax15(iT).exptCond.logPer = logPer ; %
  Data.FaulJax15(iT).exptCond.f = f ;
  Data.FaulJax15(iT).exptCond.logf = log10(f) ;

  Data.FaulJax15(iT).Results.log10_Qinv = data_Qinv(:,2) ;
  Data.FaulJax15(iT).Results.G = data_G(:,2) ;

end


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

save('ExptData.mat', 'Data')
