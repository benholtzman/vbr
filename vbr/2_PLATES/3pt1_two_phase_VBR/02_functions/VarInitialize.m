%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loads/calculates initial conditions and LABInfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vark,InitVals,LABInfo,settings,keepgoing]=...
                                  VarInitialize(Info,settings,BCs,z,dz)

%   extract some settings
    phimin=settings.phimin;

%   variables treated as constant (for now?)    
    Vark.grain = addghosts(stag(Info.init.grain)); % [m]
    Vark.grains = stag(Vark.grain); % [m]
    Vark.dg_um=Vark.grain*1e6; % [um]
    Vark.sig_MPa = addghosts(stag(Info.init.sig_MPa));
    Vark.comp=  addghosts(stag(Info.init.comp)); 
    Vark.Gu_0_GPa = addghosts(stag(Info.init.Gu_0_GPa));    
    
%   thermal properties at reference state  
    Cp_0 = addghosts(stag(Info.init.Cp_0)); % specific heat [J/kg/K]
    Rho_0 = addghosts(stag(Info.init.Rho_0)); % density [kg/m3]
    Kc_0 = addghosts(stag(Info.init.Kc_0)); % conductiviy [W/m/K]
    InitVals.Cp_0=Cp_0;InitVals.Rho_0=Rho_0;InitVals.Kc_0=Kc_0;
  
%   initial temperature
    Vark.T = addghosts(stag(Info.init.T)); % [C]
    Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);  

%   start with no phi
    Vark.phi = ones(size(Vark.T))*phimin; 
%   thermal properties at elevated P,T
    Vark = call_MaterialProps(Vark,Rho_0,Kc_0,Cp_0,settings,z,dz);  
        
%   initial water concentration
%     bulk
      Cs0 = addghosts(stag(Info.init.Cs0));    
%     get storage capacity for pargasite    
      [Vark.Strg,Vark.Kd]=pargasite_storage(z/1e3,0,settings,0,'load');    
%     put water into pargasite where it's stable
      Vark.Cs = Cs0; 
      Vark.Cparg= Cs0 .* (Vark.Strg >= Cs0);
      Bulk = Cs0 .* (Vark.Strg < Cs0); 
      
%     now calc initial Cf, Cs with initial phi guess               
      Vark = mix_Cf(Vark,Bulk,BCs,dz);
      
%   inialize CO2 -- should be zero for now....
    Vark.Cs_CO2 = addghosts(stag(Info.init.CCO2));
    Vark.Cs_CO2 = BC_setghosts(Vark.Cs_CO2,BCs.val_Cs_CO2,BCs.type_Cs_CO2,dz);
    Vark.Cf_CO2 = Vark.Cs_CO2./settings.kd_CO2;
    Vark.Cf_ORIGINAL = Vark.Cf;
    
%   Calculate solidus at initial Cf guess
    Solidus = SoLiquidus([Vark.P(2); Vark.P(2:end)],Vark.Cf,Vark.Cf_CO2,'katz');
    Vark.Tsol=Solidus.Tsol;
    Vark.Gdot = addghosts(stag(Info.init.Gdot)); % melting rate
      
      
%   get initial LAB    
%   initial viscosity   
    Vark.eta = get_VBR_visc(Vark);     
%   LAB location
    LABInfo.zLAB=0; % initialize the structure to avoid warning...
    LABInfo = find_LAB(Vark,z,settings,LABInfo);
%   add on "plume" excess
    dTex=settings.Tpot_excess-settings.Tpot; 
    zLAB=LABInfo.zLABeta; 
    dTex = dTex * (1+erf((z-zLAB)/(settings.DBL_m*2)))/2; 
    Vark.T = Vark.T + dTex; 
    Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);  
        
%   thermal properties at elevated P,T
    Vark = call_MaterialProps(Vark,Rho_0,Kc_0,Cp_0,settings,z,dz);
    InitVals.T0=Vark.T;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   volatiles and melt!  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   initializes Cf even if only doing thermal calculation so that 
%   the solidus is still calculated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
% see if max puts us above the solidus anywhere    

%   actual solidus calc (and Gdot initialization)    
    Solidus = SoLiquidus([Vark.P(2); Vark.P(2:end)],Vark.Cf,Vark.Cf_CO2,'katz');
    Vark.Tsol=Solidus.Tsol;
    Vark.Gdot = addghosts(stag(Info.init.Gdot)); % melting rate

    keepgoing = 1; % exit for time step loop
    if strcmp(settings.Flags.problem,'One_Phase_T')==0
        Molten=sum(Vark.T>Vark.Tsol);
        if Molten<=3
            disp('T < Tsol everywhere! Cannot calculate melt migration')
            disp('switching to single phase T evolution...')
            keepgoing = 1; % exit for time step loop  
            settings.Flags.problem='One_Phase_T';
        end                
    end
  
Vark.phi = ones(size(Vark.T)) * settings.phi_init; 
Vark.phi(Vark.T<Vark.Tsol)=phimin; 
Vark.phis = stag_phi(Vark,settings);
   
% Recalculate viscosity
  Vark.eta = get_VBR_visc(Vark); 

% initialize the LABInfo
  LABInfo = find_LAB(Vark,z,settings,LABInfo);
%  for time stepping:
   LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid; 
   LABInfo.zLAB1=LABInfo.zLAB; LABInfo.zLABid1=LABInfo.zLABid; 
   LABInfo.phiLAB1=LABInfo.phiLAB;LABInfo.phiLAB0=LABInfo.phiLAB;

% recalc ref muso for two-phase
  if strcmp(settings.Flags.problem,'One_Phase_T')==0
      Vark0=Vark;
      Vark0.phi=(Vark.phi>phimin)*1e-4; % at phi-crit
      eta_o = get_VBR_visc(Vark0);
      settings.muso = mean(eta_o(LABInfo.zLABid:LABInfo.zMOid));
      disp(['calculated muso: ' num2str(settings.muso)])
  end
   
  Vark.PhiFreeze = zeros(size(Vark.phi)); 
end

function zeroval=sol_Cf_finder(Cf,P,T)
Solidus = SoLiquidus(P,Cf,0,'katz');
zeroval = T - Solidus.Tsol; 
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % loads/calculates initial conditions and LABInfo
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Vark,InitVals,LABInfo,settings,keepgoing]=...
%                                   VarInitialize(Info,settings,BCs,z,dz)
% 
% %   extract some settings
%     phimin=settings.phimin;
% 
% %   variables treated as constant (for now?)    
%     Vark.grain = addghosts(stag(Info.init.grain)); % [m]
%     Vark.grains = stag(Vark.grain); % [m]
%     Vark.dg_um=Vark.grain*1e6; % [um]
%     Vark.sig_MPa = addghosts(stag(Info.init.sig_MPa));
%     Vark.comp=  addghosts(stag(Info.init.comp)); 
%     Vark.Gu_0_GPa = addghosts(stag(Info.init.Gu_0_GPa));    
% 
% %   temperature
%     Vark.T = addghosts(stag(Info.init.T)); % [C]
%     Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);  
% 
% %   thermal properties at reference state  
%     Cp_0 = addghosts(stag(Info.init.Cp_0)); % specific heat [J/kg/K]
%     Rho_0 = addghosts(stag(Info.init.Rho_0)); % density [kg/m3]
%     Kc_0 = addghosts(stag(Info.init.Kc_0)); % conductiviy [W/m/K]
%     InitVals.Cp_0=Cp_0;InitVals.Rho_0=Rho_0;InitVals.Kc_0=Kc_0;
%     InitVals.T0=Vark.T;
% 
% %   thermal properties at elevated P,T
%     Vark = call_MaterialProps(Vark,Rho_0,Kc_0,Cp_0,settings,z,dz);
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%    
% %   volatiles and melt!  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   initializes Cf even if only doing thermal calculation so that 
% %   the solidus is still calculated
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %   initial Solid water concentration
%     Cs0 = addghosts(stag(Info.init.Cs0));
% 
% %   get storage capacity for pargasite    
%     [Vark.Strg,Vark.Kd]=pargasite_storage(z/1e3,0,settings,0,'load');
% 
% %   start with melt everywhere (adjusted after solidus calc)
%     Vark.phi = addghosts(stag(Info.init.phi));
%     Vark.phi = BC_setghosts(Vark.phi,BCs.val_phi,BCs.type_phi,dz);
% 
% %   adjust solid water concentration so that (1) bulk water concentration
% %   is initially constant and (2) pargasite takes all the water (where it
% %   is stable).
% 
% %   adjust Cs with mixing of incompatible element:
%     Vark.Cs= Cs0./((1./Vark.Kd).*Vark.phi + (1-Vark.phi));
%     Vark.Cs = BC_setghosts(Vark.Cs,BCs.val_Cs,BCs.type_Cs,dz);
% 
% %   put water into pargasite where it's stable
%     Vark.Cparg= Cs0 .* (Vark.Strg > Vark.Cs); 
%     Vark.Cs(Vark.Strg > Vark.Cs)=0;
% 
% %   calculate Cf and Bulk
%     Vark.Cf = Vark.Cs./Vark.Kd; % this would be max Cf... prob should calc initial phi s.t. Tsol(Cf)=T with Cf<Cs/kd
%     Vark.Cbulk= Vark.Cf.*Vark.phi + Vark.Cs .* (1-Vark.phi);
%     
% %   (should do the same for CO2, but this isn't used right now...)    
%     Vark.Cs_CO2 = addghosts(stag(Info.init.CCO2));
%     Vark.Cs_CO2 = BC_setghosts(Vark.Cs_CO2,BCs.val_Cs_CO2,BCs.type_Cs_CO2,dz);
%     Vark.Cf_CO2 = Vark.Cs_CO2./settings.kd_CO2;
%     Vark.Cf_ORIGINAL = Vark.Cf;
% 
% % now calculate solidus and adjust phi appropriately
% %   solidus (calculate even if not solving for phi evolution)
%     Solidus = SoLiquidus([Vark.P(2); Vark.P(2:end)],Vark.Cf,Vark.Cf_CO2,'katz');
%     Vark.Tsol=Solidus.Tsol;
%     Vark.Gdot = addghosts(stag(Info.init.Gdot)); % melting rate
% 
% %   adjust phi --> phimin where below solidus
%     Vark.phi(Vark.T<Vark.Tsol) = phimin;          
%     Vark.phis = stag_phi(Vark,settings);
%    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % add on excess temperature below LAB
% % then re-calculate initial phi and Cf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % initial viscosity   
%   Vark.eta = get_VBR_visc(Vark); 
%     
% % LAB location
%   LABInfo.zLAB=0; % initialize the structure to avoid warning...
%   LABInfo = find_LAB(Vark,z,settings,LABInfo);
% 
% % add on "plume" excess
%   dTex=settings.Tpot_excess-settings.Tpot; 
%   zLAB=LABInfo.zLABeta; 
%   dTex = dTex * (1+erf((z-zLAB)/(settings.DBL_m*2)))/2; 
%   Vark.T = Vark.T + dTex; 
%   Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);  
% 
% % recalculate density, conductivity, etc.  
%   Vark = call_MaterialProps(Vark,Rho_0,Kc_0,Cp_0,settings,z,dz);
% 
% % Recalculate solidus, phi and viscosity  
%   Solidus = SoLiquidus([Vark.P(2); Vark.P(2:end)],Vark.Cf,Vark.Cf_CO2,'katz');
%   Vark.Tsol=Solidus.Tsol;
%   Vark.phi = settings.phi_init*ones(size(Vark.phi));
%   Vark.phi = BC_setghosts(Vark.phi,BCs.val_phi,BCs.type_phi,dz);
%   Vark.phis = stag_phi(Vark,settings);
%   Vark.phi(Vark.T<Vark.Tsol) = phimin;          
%   Vark.eta = get_VBR_visc(Vark); 
%  
% % Recalculate water content
%   Vark.Cs= Cs0./((1./Vark.Kd).*Vark.phi + (1-Vark.phi));
%   Vark.Cparg= Cs0 .* (Vark.Strg > Vark.Cs); 
%   Vark.Cs(Vark.Strg > Vark.Cs)=0;
%   Vark.Cs(Vark.Strg > Vark.Cparg)=0;
%   Vark.Cs = BC_setghosts(Vark.Cs,BCs.val_Cs,BCs.type_Cs,dz);
%   Vark.Cf = Vark.Cs./Vark.Kd;
%   Vark.Cbulk= Vark.Cf.*Vark.phi + Vark.Cs .* (1-Vark.phi);
%     plot(Vark.Cf); hold all; 
% % initialize the LABInfo
%   LABInfo = find_LAB(Vark,z,settings,LABInfo);
% %  for time stepping:
%    LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid; 
%    LABInfo.zLAB1=LABInfo.zLAB; LABInfo.zLABid1=LABInfo.zLABid; 
%    LABInfo.phiLAB1=LABInfo.phiLAB;LABInfo.phiLAB0=LABInfo.phiLAB;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % check if there is melt %
% %%%%%%%%%%%%%%%%%%%%%%%%%%
%   keepgoing = 1; % exit for time step loop
%   if strcmp(settings.Flags.problem,'One_Phase_T')==0
%      Molten=sum(Vark.T>Vark.Tsol); 
%      if Molten==0
%         disp('T < Tsol everywhere! Cannot calculate melt migration')
%         disp('Finishing initialization then exiting...')	
%         keepgoing = 0; % exit for time step loop
%      end
%   end
% 
% end
