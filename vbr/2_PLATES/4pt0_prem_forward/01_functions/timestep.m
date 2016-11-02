%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vark1,resid,tnow_s_1,LABInfo] = ...
       timestep(Vark,tnow_s,LABInfo,settings,z,dz,IVals,BCs)
                  
   dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]   
   TempUpdate= 1-strcmp(settings.Flags.TempUpdate,'DiffusionOnly');
   cfl = settings.cfl; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Calculate temperature fluxes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  adiabatic heat flux
   Fluxes.adi = Vark.Vbgz.* settings.dTdz_ad; 
   
%  diffusive heat flux   
%    A_T = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,Vark.rho,Vark.cp,Vark.Kc);    
%    dTdt.diff(:) = (A_T(:,:)) * Vark.T(:); % diffusive
   Fluxes.diff = make_dTdt_diff(z,Vark.rho,Vark.cp,Vark.Kc,Vark.T);% * TempUpdate; % CJH
   Fluxes.diff(1) = 0; Fluxes.diff(end) = 0; 
    
%  calculate stable thermal step
   dt_TDiff = 0.5*dz*dz/max(Vark.Kc./Vark.rho./Vark.cp);   
   dt_TAdv = dtCFL(Vark.Vbgzs,dz,cfl,dt_max);
   dt_adi = 0.1 * min(abs(Vark.T./Fluxes.adi));
   dt=min([dt_TDiff dt_TAdv dt_adi]);
   
%  advective heat fluxes  
   Fluxes.adv = advection_driver(Vark.T,Vark.Vbgzs,dz,dt,'VL_nc');
   
% calculate full operator
   Fluxes.Total =  (Fluxes.diff + Fluxes.adi*TempUpdate + Fluxes.adv*TempUpdate);

%%%%%%%%%%%%%%%%%%%%%%%
%% update temperature %
%%%%%%%%%%%%%%%%%%%%%%%
     
%   save old values
    Old.T = Vark.T;
    
%   update  
    nz = numel(z); 
    Vark.T(2:nz-1) = Vark.T(2:nz-1) + Fluxes.Total(2:nz-1)*dt; 
    
%   set BCs    
%     dTdz = calc_dTdz_LAB(Vark,settings); 
%     [BCs]=init_BCs(BCs,'T','zmax','neumann',dTdz);
    Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);      
    tnow_s_1 = tnow_s + dt;
    
%   calculate residuals      
    resid.T = max(abs(Vark.T(2:nz-1) - Old.T(2:nz-1))./(Old.T(2:nz-1)));        
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% calculate "static" quantities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate solidus, set phi to phi_min above solidus
  Solidus = SoLiquidus([Vark.P(2); Vark.P(2:end)],Vark.Cf_H2O,Vark.Cf_CO2,'katz');
  Vark.Tsol=Solidus.Tsol;
  Vark.phi = ones(size(Vark.T)) * settings.phi_min; 
  Vark.phi(Vark.T<Vark.Tsol)=0;
  
% density, thermal properties: 
  Vark = call_MaterialProps(Vark,IVals.Rho_0,IVals.Kc_0,IVals.Cp_0,...
                            settings,z,dz);  
  Vark.eta = get_VBR_visc(Vark); 
  
% current LAB location
  [LABInfo] = find_LAB(Vark,z,settings,LABInfo);    

% output the final arrays
  Vark1 = Vark; 
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dt]=dtCFL(vf,dz,cfl,dt_max)
   % calculate time step      
     vfmax = max(abs(vf)); 
     dt = cfl*dz/(vfmax+1e-20); % porous flow max time step
     mindt = dt_max *1e6 * 365 * 24 * 3600; 
     dt = dt*(vfmax~=0)+mindt*(vfmax==0);  % control for vf == 0     
end

