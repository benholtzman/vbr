%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vark1,resid,tnow_s_1,dt,dt_structs] = timestep(Vark,tnow_s,zLAB,...
                      phiLAB,zLABid,settings,nz,z,zs,dz,IVals,BCs,cfl,in_i)
   dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]
   phimin =  settings.phimin; 
   phimax =  settings.phimax; 
   Vbgs = -settings.Vbg / (1e2 * 365 * 24 * 3600); % [m/s] for mass conservation
   DBL_m = settings.DBL_m; % mechanical boundary layer thickness at LAB
   VbgFlag = settings.Flags.VbzFlag; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate fluxes and source/sink terms %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dTdt  structure with each component:
%       .adi = adiabatic heat flux
%       .diff = diffusive heat flux 
%       .Gdot = latent heat source/sink
%       .adv  = advective heat flux
%    dphidt  conservative melt flux  
%    Gdot    melting rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  zero everything 
   [dTdt,dphidt]=initialize_dT_dphi(nz);
   
%  non-advective fluxes      
   dTdt.adi(in_i) = Vark.Tot(in_i).* settings.dTdz_ad; % adiabatic
   A_T = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,Vark.rho,Vark.cp,Vark.Kc); 
   dTdt.diff(:) = (A_T(:,:)) * Vark.T(:); % diffusive
   dTdt.diff(1) = 0; dTdt.diff(end) = 0; 
   dCfdt = dCfdt_init(Vark,settings,dz); 
      
%  calculate stable time steps, use minimum
   dt_TDiff = 0.5*dz*dz/max(Vark.Kc./Vark.rho./Vark.cp);   
   dt_TAdv = dtCFL(Vark.Tots,dz,cfl,dt_max);
   dt_CfDiff = 0.5*dz*dz/dCfdt.maxdiff;
   dt_CfCFL=dtCFL(dCfdt.Veff,dz,cfl,dt_max);
   dt_phiCFL=dtCFL(Vark.vfs,dz,cfl,dt_max);
%    dt=min([dt_CfDiff dt_CfCFL dt_phiCFL dt_TDiff dt_TAdv ]);
   dt=min([dt_TDiff dt_TAdv ]);
% calculate advective fluxes
  [dTdt,dphidt,dCfdt]=advective_fluxes(Vark,dz,dt,dCfdt,dphidt,dTdt);

% calculate melting/crystallization rate
   Vgdot = Vark.Vbgz; %Vark.Vs;
   Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLAB,dphidt);   
   [dphidt,dTdt,dCfdt]=SourceTerms(Vark,settings,dCfdt,dphidt,dTdt);

% limit how much phi can change in a step
   fstep = 0.02;
   phiproject=abs((dphidt.adv(in_i)*dt + Vark.Gdot(in_i)*dt));
   if max(abs(phiproject)) > fstep
     dt_m = min(abs((fstep)./(dphidt.adv(in_i) + Vark.Gdot(in_i))));
     dt = min([dt_m dt]); 
%    recalculate the terms that depend on time step length
     [dTdt,dphidt,dCfdt]=advective_fluxes(Vark,dz,dt,dCfdt,dphidt,dTdt);
     Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLAB,dphidt);
     [dphidt,dTdt,dCfdt]=SourceTerms(Vark,settings,dCfdt,dphidt,dTdt);
   end
   
% calculate full operators
   dTdt.Total =  (dTdt.diff + dTdt.adi + dTdt.adv + 0*dTdt.Gdot); 
%    dCfdt.diff = dCfdt.diff*0; dCfdt.adv = 0; dCfdt.Gdot = 0; %CJH
   dCfdt.Total = 0*(dCfdt.diff + dCfdt.adv + dCfdt.Gdot);
%    dphidt.Total = (dphidt.adv.*(Vark.T > Solidus.Tsol)  + dphidt.Gdot); 
   dphidt.Total = 0*(dphidt.adv  + dphidt.Gdot + 0*Vark.Sd_Gdot); 
   dt_structs.dphidt = dphidt;
   dt_structs.dCfdt = dCfdt;
   dt_structs.dTdt = dTdt;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update the time dependent variables %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
%   save old values
    Old.T = Vark.T; Old.phi = Vark.phi; Old.Cf=Vark.Cf; 
    
%   update     
    Vark.T(in_i) = Vark.T(in_i) + dTdt.Total(in_i)*dt;  
    Vark.phi(in_i) = Vark.phi(in_i) + dphidt.Total(in_i)*dt;
    Vark.Cf(in_i) = Vark.Cf(in_i) + dCfdt.Total(in_i)*dt;
    tnow_s_1 = tnow_s + dt;
    
%   set BCs and phi-limits  
%     T
      Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);    
%     phi  
      Vark.phi(Vark.phi<phimin)=phimin;
      Vark.phi(Vark.phi>phimax)=phimax;      
      Vark.phi = BC_setghosts(Vark.phi,BCs.val_phi,BCs.type_phi,dz);                   
%     Cf (and Cs)               
      Vark.Cs = Vark.Cf.*settings.kd_H2O;
      Vark.Cs = BC_setghosts(Vark.Cs,BCs.val_Cs,BCs.type_Cs,dz);
      Vark.Cf = BC_setghosts(Vark.Cf,BCs.val_Cs./settings.kd_H2O,BCs.type_Cs,dz);

%  calculate residuals      
    resid.T = max(abs(Vark.T(in_i) - Old.T(in_i))./(Old.T(in_i)));
    resid.phi = max(abs(Old.phi(in_i) - Vark.phi(in_i))./Old.phi(in_i));
    resid.Cf = max(abs(Old.Cf(in_i) - Vark.Cf(in_i))./Old.Cf(in_i));
    resid.maxresid = max([resid.phi resid.T]);       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% calculate "static" quantities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% density, thermal properties: 
  Vark = call_MaterialProps(Vark,IVals.Rho_0,IVals.Kc_0,IVals.Cp_0,settings,z,dz);  
   
% melt system:   
%    calculate pemeability        
     Vark.phis = stag_phi(Vark,settings);
     Vark=permeability(Vark,settings);
%    calculate surface tension coefficient on staggered grid
     SurfCoff=SurfaceTension(Vark,settings); 
%    update LAB location
     LABInfo = find_LAB(Vark,z,settings,LABInfo);

     if strcmp(VbgFlag,'variable')
         Vark.Vbgzs=upwelling_velocity(Vbgs,zs,zLAB,DBL_m);
         Vark.Vbgz = addghosts(stag(Vark.Vbgzs));
     elseif strcmp(VbgFlag,'con_asthen')
         Vark.Vbgzs = Vbgs.*ones(size(Vark.phis));
         Vark.Vbgzs(zs<zLAB)=0;
         Vark.Vbgz = addghosts(stag(Vark.Vbgzs));
     end
     Vark = calc_Sd(Vark,settings,z,zs,zLAB,Vark.Pex);           
     Vark.Sd_Gdot = calc_Sd_sink(Vark.Sds,z,dz,zLABid); 
     
      if strcmp(settings.Flags.problem,'Two_Phase_Sd')
%        extract only the section below the solidus-liquidus intersection
         [Molten,phiLAB] = prep_molten(Vark,SurfCoff,Vark.Vbgzs,zLABid,settings);          
%        calculate the boundary conditions!          
         Sd_type = BCs.type_P; Sd_val = BCs.val_P; 
         Sd_type(1) = 2; % it's a flux!
         [Sd_val(1),Sdflux] = calc_Sd_dPdz(Molten.phi(2),Molten.perms(1),...
                     Molten.SurfTen(2),settings,dz);          
%        and now do the calculation         
         [Molten.Pex,A,src] = P_solve_gen(Molten.phi,Molten.phis,Molten.perms,...
                      Molten.mus,Molten.SurfTen,0,dz,settings,Sd_val,Sd_type);
%        Get flux and velocities    
         Molten = darcy(Molten,Molten.Vbgs,Molten.SurfTen,0,dz,settings);          
%        store in global array
         Vark = store_molten(Vark,Molten,zLABid);         
     else
         [Vark.Pex,A,src] = P_solve_gen(Vark.phi,Vark.phis,Vark.perms,...
                                Vark.mus,SurfCoff,0,dz,settings,BCs.val_P,BCs.type_P);        
    %    Get flux and velocities    
         Vark = darcy(Vark,Vark.Vbgzs,SurfCoff,0,dz,settings);      
      end
          
      Vark1 = Vark;       
end
