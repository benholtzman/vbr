function PLOT_bestfit_PTdep_Hspace
close all;

% load Box data
  Box_dir = '../../../../0_BOXES/';
  
  
  
% store and save the best fit indeces
 
fsize = 12;
ylimits = [0 450];
 
 
 for iBox = 1:3
    if iBox == 1
        Box_name ='Box_halfspace_con_ad4_150225_VBR.mat';
        Fig1 = figure('color',[1 1 1]);
        clr = 'k';
        dnm = 'Constant';
    elseif iBox == 2
        Box_name ='Box_halfspace_P_dep_ad4_150225_VBR.mat';
        clr = 'b';
        dnm = 'f(P)';
    elseif iBox == 3
        Box_name ='Box_halfspace_PT_dep_ad4_150225_VBR.mat';
        clr = 'r';
        dnm = 'f(P,T)';
    end
   % load it     
     Box_name = [Box_dir Box_name]
     load(Box_name)
   % extract data
     iV1 = BestBox.iVar1;
     iV2 = BestBox.iVar2;
     tsnap= BestBox.time_index;
     
     Z_km = Box(iV1,iV2).Movie.info.Z_km;
     Rho = Box(iV1,iV2).Movie.Frames(tsnap).rho;
     Kc = Box(iV1,iV2).Movie.Frames(tsnap).Kc;
     
     if iBox == 3 % (Cp wasn't being save for these runs... fixed that, but
                  % need to recalculate it here for now)
         T = Box(iV1,iV2).Movie.Frames(tsnap).T;
         Cp = SpecificHeat(T,0.9);
     else         
         cp_olv = 1100;		% Specific heat J /(kg K) PER UNIT MASS
         cp_crust = 800;		% (????) Specific heat J /(kg K) (only used if Kc_rho_Cp_type = con)
         Cp = Rho./Rho*cp_olv; 
         Cp(Z_km<=Box(iV1,iV2).Movie.info.Z_moho)=cp_crust;
     end
     
     
     subplot(1,4,1)
         if iBox > 1; hold on; end
         plot(Rho,Z_km,clr,'displayname',dnm)
         if iBox > 1; hold off; end
         xlabel('Density (kg m^{-3})','fontsize',fsize); 
         ylabel('Depth (km)','fontsize',fsize)
         set(gca,'Ydir','rev')
         
     subplot(1,4,2)
         if iBox > 1; hold on; end
         plot(Kc,Z_km,clr,'displayname',dnm)
         if iBox > 1; hold off; end
         xlabel('Conductivity (W m^{-1} K^{-1})','fontsize',fsize); 
         ylabel('Depth (km)','fontsize',fsize)
         set(gca,'Ydir','rev')  
 ylim(ylimits)
         
      subplot(1,4,3)
         if iBox > 1; hold on; end
         plot(Cp,Z_km,clr,'displayname',dnm)
         if iBox > 1; hold off; end
         xlabel('Specific Heat (J kg^{-1} K^{-1})','fontsize',fsize); 
         ylabel('Depth (km)','fontsize',fsize)
         set(gca,'Ydir','rev')       
 ylim(ylimits)
         
      subplot(1,4,4)
         if iBox > 1; hold on; end
         semilogx(Kc./Cp./Rho,Z_km,clr,'displayname',dnm)
         if iBox > 1; hold off; end
         xlim([1e-7 1e-6])
         xlabel('Diffusivity (m^2 s^{-1})','fontsize',fsize); 
         ylabel('Depth (km)','fontsize',fsize)
         set(gca,'Ydir','rev')     
 ylim(ylimits)
     
     
 end
 
 
 load ../../PREM/PREM.mat
 subplot(1,4,1)
 hold on
 % take off water/sed, adjust depth to 0
 PREM_rho = Density(Density>2000);
 PREM_z = (max(Radius)-Radius)/1e3;
 PREM_z = PREM_z(Density>2000);
 PREM_z = PREM_z - min(PREM_z);

 plot(PREM_rho,PREM_z,'--k','displayname','PREM')
 ylim(ylimits)
 xlim([2 4]*1e3)
 
end

function Cp = SpecificHeat(T,FracFo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates specific heat as a function of temperature using a polynomial.
% 
% Input
%  
%   T       temperature [K], array or scalar
%   FracFo  volume fraction forsterite
%
% Output
%
%   Cp      Specific Heat [J kg^-1 K^-1]
%
% Papers report heat capacity (J/mol/K) as polynomial functions of
% temperature with coefficients determined by fitting. To get to specific
% heat (J/kg/K), I divide the reported coefficients by molecular weight of 
% the minerals. The function form is typically:
% 
%  Cp = B(1) + B(2)*T + B(3)*T^2 + B(4)*T^-2 + B(5)*T^-3 + ...
%       B(6)*T^-0.5 + B(7)/T
%
% In this implementation, the array B initially has two rows, with values 
% for forsterite (Fo) in the first row and fayalite (Fa) in the second row. 
% The two are then linearly weighted by the fraction of forsterite in the 
% mantle before calculating Cp. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need molecular weights of forsterite (Mg2Si04) and fayalite (Fe2SiO4): 
% (reminder, 1 atomic unit = 1 g / mol = 1e-3 kg / mol)
  wt_Fo = (24.305*2+28.085 + 15.999*4)/1000; % [kg/mol]
  wt_Fa = (55.845*2+28.085 + 15.999*4)/1000; % [kg/mol]  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE WHICH FIT TO USE %
% (pretty much the same)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% Fitting 1: Saxena, S. K. "Earth mineralogical model: Gibbs free energy %
% minimization computation in the system MgO-FeO-SiO2." Geochimica et    %
% cosmochimica acta 60.13 (1996): 2379-2395.                             %
%                                                                        %
% Fitting 2: Berman, R. G., and L. Ya Aranovich. "Optimized standard     %
% state and solution properties of minerals." Contributions to Mineralogy% 
% and Petrology 126.1-2 (1996): 1-24.                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % FITTING 1 
%    B = [165.80 0.1855*1e-1 0.0 -3971000 0.28161*1e9 0.0 -5610.0; ...
%         104.45 0.4032*1e-2 0.0 2947000 -.3941*1e9 0.0 -173330.0];
   
% FITTING 2: Berman and Aranovich, Contrib. Mineral. Petrol. 1996
   B = zeros(2,7); 
   B(1,1)=233.18; B(1,6)=-1801.6;B(1,5)=-26.794*1e7; 
   B(2,1)=252;B(2,6)=-2013.7; B(2,5)=-6.219*1e7;         

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% convert from heat capacity to specific heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   B(1,:) = B(1,:) / wt_Fo; % [J/mol/K to J/kg/K]
   B(2,:) = B(2,:) / wt_Fa; % [J/mol/K to J/kg/K]
   B = B(1,:)*FracFo + B(2,:)*(1-FracFo); 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% calculate heat capacity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Cp = B(1) + B(2)*T + B(3)*T.^2 + B(4)*T.^(-2) + B(5)*T.^(-3) + ...
       B(6)*T.^(-0.5) + B(7)./T;
   
end
