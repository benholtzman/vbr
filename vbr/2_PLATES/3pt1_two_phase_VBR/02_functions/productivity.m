%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads melt productivity curves from isentropic melting. 
%
% input
%    zNum   depth 
%    CompID composition ID: 1 for peridotite batch, 2 for peridotite
%                           fractional, 3 for peridotite-pyroxenite
%    Frac   if CompID = 1 or 2, Frac is the ppm H2O in source
%           if CompID = 3, Frac is the percent of pyroxenite component in
%           source
%    Tpot   mantle potential temperature [C]
% output
%    dFdP   melt productivity interpolated to zNum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFdP = productivity(zNum,CompID,Frac,Tpot)

if CompID == 1 % peridotite batch melting (works for dry)
    lnm=['../../dFdP/data_peridotite/Perid_Batch_H2O_' num2str(Frac) '_Tpot_' num2str(Tpot) '.mat'];            
elseif CompID == 2 % peridotite pure fractional melting (works for dry)
    lnm=['../../dFdP/data_peridotite/Perid_Frac_H2O_' num2str(Frac) '_Tpot_' num2str(Tpot) '.mat'];            
elseif CompID == 3 % anhydrous peridotite-pyroxenite melting 
    lnm=['../../dFdP/data_pyrxper/Mixture_PhiPyx_' num2str(Frac) '_Tpot_' num2str(Tpot) '.mat'];   
end

if (exist(lnm,'file')==0) 
      display('dFdP not calculated, doing so now')
      addpath('../../dFdP/')
      
%   set variables used by all dFdP drivers
      svenme = lnm;
      InVars.Mcpx = .25; % CPX weight fraction (not really used...)
      InVars.Zstart = 300; % depth to start at [km]
      InVars.Zstop = 30; % depth to stop at [km]
      InVars.nz =600; % number of z points
      InVars.Tpot = Tpot; % mantle potential temperature
      % InVars.To = 0.6*InVars.Zstart+InVars.Tpot; % initial temperature at depth [oC]
%       InVars.To = InVars.Tpot*exp(InVars.Zstart*1e3*9.8*3300*30*1e-6/3300/1e3); % initial temperature at depth [oC]
      InVars.RKOrder=4; % 2 or 4 for 2nd or 4th order runge-kutta integration
% things for adiabat
%     X = [pyroxenite from shorttle; peridotite from katz]
      Cp = [1140; 1000];   
      al_s = [40*1e-6; 40*1e-6];
      rho_s = [3300; 3300];
    
%   set driver specific settings and call driver
      if CompID == 1
          H2OPPm = Frac; % bulk water, ppm.
          InVars.H2OBulk = H2OPPm/1e6*100; % bulk water wt%
          
%        average expansion and specific heat coefficients
          phi = [0; 1];
          expan = sum(phi.*al_s./rho_s); 
          Cpbar = sum(phi.*Cp);      
          coff = expan/Cpbar; 
          InVars.To = Tpot*exp(InVars.Zstart*3300*9.8*1e3*coff); % adiabatic geotherm extrapolation 
          
          [Vars]=peridotite_melting(InVars);
      elseif CompID == 2
          H2OPPm = Frac;
          InVars.H2OBulk = H2OPPm/1e6*100; % bulk water wt%
          InVars.dFextract = .0000001;
          
%        average expansion and specific heat coefficients
          phi = [0; 1];
          expan = sum(phi.*al_s./rho_s); 
          Cpbar = sum(phi.*Cp);      
          coff = expan/Cpbar; 
          InVars.To = Tpot*exp(InVars.Zstart*3300*9.8*1e3*coff); % adiabatic geotherm extrapolation
                    
          [Vars]=peridotite_fractional_melting(InVars);
      elseif CompID == 3
          Pyx = Frac/100; 
          InVars.PhiPyx = Pyx*(Pyx>0)*(Pyx<1)+0.00001*(Pyx==0)+0.999*(Pyx==1);
          InVars.PhiPer = 1-InVars.PhiPyx;
          InVars.H2OBulk = 0.0; % bulk water wt%
          
%        average expansion and specific heat coefficients
          phi = [InVars.PhiPyx; InVars.PhiPer];
          expan = sum(phi.*al_s./rho_s); 
          Cpbar = sum(phi.*Cp);      
          coff = expan/Cpbar; 
          InVars.To = Tpot*exp(InVars.Zstart*3300*9.8*1e3*coff); % adiabatic geotherm extrapolation 
          
%          [Vars]=mixed_melting(InVars);
          [Vars]=two_component_melting(InVars);
      end
      
%   save it so we don't have to do this every time
      save(svenme,'Vars','InVars');
      rmpath('../../dFdP')
      display(['dFdP calculated for' num2str([CompID Frac Tpot])])
else
    load(lnm)
end
    
  dFdPcalcd=Vars.dFdP(1*(CompID<3)+3*(CompID==3),:); % use the aggregate melt productivity if mixture
  zdFdP=Vars.z;
% actually do what we came for:
  dFdP=interp1(zdFdP,dFdPcalcd,zNum); % in 1 / GPa
end


