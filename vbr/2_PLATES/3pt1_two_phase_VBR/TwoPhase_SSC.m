%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      
%%                        TwoPhase.m 
%%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1D finite-volume solution to two-phase equations for conservation of
%% mass, momentum and energy. Solution is dimensional.  
%% 
%% Assumptions:
%%  - uniform mesh
%%  - no pressure gradients from solid shear
%%  - includes surface tension term (not yet...)
%%  - melting rate = material derivative of F(P,T,Tsol,Composition)
%%
%% Many initial conditions and material settings are controlled by the 
%% functions in 01_functions/ and 01_init_conds_and_material. 
%% These functions are all called by driver.m before calling TwoPhase. 
%%
%% For finite volume, scalar quantities (pressure, melt fraction) are 
%% solved at cell centers while velocities are defined on cell edges.
%%
%% Initial conditions are given as input on cell edges. 
%% 
%% Final Output is output on cell edges in 2D arrays where each column 
%% corresponds to a single timestep. i.e., phi(:,3) is the melt fraction
%% profile at time step 3. Variables are stored in the Vars structure:
%% 
%% Vars.
%%    P        solid pressure [Pa]
%%    Pex      excess pressure [Pa]
%%    phi      porosity
%%    S        darcy flux [m s^-1]
%%    V        solid velocities [m s^-1]
%%    vf       fluid velocities [m s^-1]
%%    T        temperature [C]
%%    Gdot     melting rate [s^-1]
%%    cp       specific heat capacity [J kg^-1 K^-1]
%%    rho      density (solid?) [kg m^-3]
%%    Kc       thermal conductivity [W m^-2]
%%    dg_um    grain size [micrometers]
%%    sig_MPa  stress [MPa] 
%%    eta      shear viscosity [Pa s] (calculated via VBR, max 1e26)
%%    Cs       solid tracer concentration (water in this case) [PPM]
%%    Cf       fluid tracer concentration [PPM]
%%    comp     compositional weighting function [crust = 0, mantle = 1]
%%    Gu_0_GPa unrelaxed elastic shear moduli at reference P,T [GPa]
%%
%%    Need to update the list! Parg!
%%
%% Info.
%%    z        depth array [m]
%%    z_km     depth array [km]
%%    t        array with times corresponding to the columns in
%%             the solution arrays [s]
%%    t_Myr    same as t except [Myr]
%%    init.    structure containing initial conditions
%%    final    end state flag (0 is clean exit, 1 is steady state)
%%    ssresid  final max resid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vars,Info]=TwoPhase_SSC(Info,settings)

%% --------------------------------------------------------------------- %%
%% ----------               Initialization                      -------- %%
%% --------------------------------------------------------------------- %%
  tinit = tic; % initialize time counter
  addpath ./02_functions % load in functions
  addpath ./03_plotting % load in functions
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% load in settings and things %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% space   
  dz = settings.Zinfo.dz_m; % mesh spacing
  zs = settings.Zinfo.z_m; % vertical coordinate of cell edges
% time
  nt = settings.nt; % time steps to go! 
  outk = settings.outk; % frequency of output (output every outk steps)
  outn = nt/outk; % number of timesteps to save    
% solution settings
  cfl = 0.4;      % advection stability criteria  
%   dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]
  ss_tol = settings.sstol; % steady state tolerance for normalized phi residual. 
  prog_pause = 0.001; 
% Initialize flags and counters
  VbgFlag = settings.Flags.VbzFlag; 
  k = 1;         % time step count
  tnow_s = 0;      % current model time
  t_max_Myrs = settings.t_max_Myrs; % max time [Myrs]
  t_s_to_Myrs = 3600*24*365*1e6; % seconds --> Myr factor
% Boundary Conditions 
  BCs = Info.BCs; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% mesh and initial values %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% input z is the staggered (cell edge) z, need cell-centered
  z = stag(zs);     
  z = [z(1)-dz; z(:); z(end)+dz];
  nz = numel(z);
  in_i = 2:nz-1; % node mask for nonghosts
  Info.z = z; 
  Info.z_km = z/1e3;

% initialize output structures:
%    Vars: output structure, Vars.variable(depth,time)
%    Info: structure with all the settings, 1D results (e.g., LAB depth)
%    kk: output index counter
  [Vars,Info,kk]=Varstruct(Info,z,outn,settings);

% initialize working structures
%    Vark: current time step values
%    InitVals: reference state for material properties 
%    LABInfo: structure with LAB and solidus-geotherm depth, indeces
%    keepgoing: time loop kill flag
  [Vark,InitVals,LABInfo,settings,keepgoing] = ...
                                 VarInitialize(Info,settings,BCs,z,dz);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Initial Dependent Fields %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate pemeability        
  Vark=permeability(Vark,settings);
% calculate surface tension coefficient on staggered grid
  SurfCoff=SurfaceTension(Vark,settings); 
% background upwelling velocity
  Vbgs=-settings.Vbg/(1e2 * 365 * 24 * 3600); % [m/s] 
  Vark.Vbgzs = Vbgs.*ones(size(Vark.phis));
  if strcmp(settings.Flags.LABdef,'visc')==1
         zLAB_Vbg=LABInfo.zLABeta;
  elseif strcmp(settings.Flags.LABdef,'phi')==1
         zLAB_Vbg=LABInfo.zLAB;
  end  
  if strcmp(VbgFlag,'variable') || strcmp(VbgFlag,'var_z_con_t')
      Vark.Vbgzs=upwelling_velocity(Vbgs,zs,zLAB_Vbg,settings.DBL_m);
  elseif strcmp(VbgFlag,'con_asthen')
      Vark.Vbgzs(zs<zLAB_Vbg)=0;
  end
  Vark.Vbgz = addghosts(stag(Vark.Vbgzs));

  if keepgoing * strcmp(settings.Flags.problem,'One_Phase_T')==0
% Solve for P
  [Vark.Pex,~,~] = P_solve_gen(Vark.phi,Vark.phis,Vark.perms,Vark.eta,...
                       SurfCoff,dz,settings,BCs.val_P,BCs.type_P); 
  else
      Vark.Pex=zeros(size(Vark.T)); 
  end
                   
% Get flux and velocities:      
  Vark = darcy(Vark,Vark.Vbgzs,SurfCoff,dz,settings); 

% run some diffusion iterations on T to smooth out initial condition
  settings.Flags.TempUpdate='DiffusionOnly'; tnow_init = 0;
%   figure
  for it_k0 = 1:50
  [Vark,resid,tnow_init,dt,LABInfo0,Extruded] = split_timestep(Vark,tnow_init,LABInfo,...
                                settings,nz,z,zs,dz,InitVals,BCs,cfl,in_i,0);
% %    plot(Vark.phi,'color',[it_k0/100 0 0]); hold on
%    subplot(2,1,1); hold on
%    plot(Vark.T,'color',[it_k0/50 0 0]); 
%    
%    subplot(2,1,2); hold on
%    plot(Vark.Gdot,'color',[it_k0/50 0 0]); 
%    pause(0.03)
  end
  
  
% save initial step   
  [Vars,Info,kk]=savevar(Vars,Vark,Info,tnow_s,0,k,kk,LABInfo);

% plot initial step  
  if strcmp(settings.Flags.progress_plot,'yes')
    Vark.fprog.fighand=[]; 
    Vark.fprog.x_width=5 ;Vark.fprog.y_width=8; % [inches
    Vark.fprog.prog_pause=prog_pause; 
    Vark.fprog = progress_plot(Vark,z,Vark.fprog);
  end
  
%% --------------------------------------------------------------------- %%
%% ---------- Solve Forward Problem (time stepping starts here) -------- %%
%% --------------------------------------------------------------------- %%
settings.Cbulk0=Vark.Cbulk;   
settings.Flags.TempUpdate='yes'; 
LABInfo.lag_steps=0; Extruded = 0; 
while keepgoing == 1 && k <= nt
    k = k + 1;
    
%%%%%%%%%%%%%%%
%% Time step %%
%%%%%%%%%%%%%%%

%   [Vark,resid,tnow_s,dt,dt_structs]= timestep(Vark,tnow_s,zLAB,phiLAB,zLABid,...
%                                 settings,nz,z,zs,dz,InitVals,BCs,cfl,in_i);
[Vark,resid,tnow_s,dt,LABInfo,Extruded] = split_timestep(Vark,tnow_s,LABInfo,...
                                settings,nz,z,zs,dz,InitVals,BCs,cfl,in_i,Extruded);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output and error/ss check %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if mod(k,outk)==0;                         
        [Vars,Info,kk]=savevar(Vars,Vark,Info,tnow_s,Extruded,k,kk,LABInfo);
        if strcmp(settings.Flags.problem,'One_Phase_T')==0
         disp(' ');disp(['     max phi: ' num2str(max(Vark.phi))])
         disp(['     phi residual: ' num2str(resid.phi)])             
        end
        if strcmp(settings.Flags.H2O,'variable')==1
         disp(['     Cf residual: ' num2str(resid.Cf)]) 
        end
        disp(['     T residual: ' num2str(resid.T)])
        disp(['     zsol [km]: ' num2str(LABInfo.zLAB/1e3)])        
        disp(['     zLAB [km]: ' num2str(LABInfo.zLABeta/1e3)])   
        disp(['     Extruded volume thickness [m]: ' num2str(Extruded)])
    end 
    if isreal(resid.T)==0 || isreal(dt) == 0
        disp('   wtf')
    end
    
    [keepgoing,Info] = check_stop_run(keepgoing,Info,resid,ss_tol,tnow_s,...
                                    t_max_Myrs,LABInfo,settings,Vark,k <= nt);
    

    Molten=sum(Vark.T>Vark.Tsol);
end
%% --------------------------------------------------------------------- %%
%% ----------               Endgame                             -------- %%
%% --------------------------------------------------------------------- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% final save and cleanup %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% final variable save
  disp(' '); disp('Calculations complete, wrapping up...')  
  [Vars,Info,kk]=savevar(Vars,Vark,Info,tnow_s,Extruded,k,kk,LABInfo);
  [Vars,Info]=varfinal(Vars,Info,kk); % removes unfilled columns  
  Info.tMyrs=Info.t/3600/24/365/1e6;
  Info.ssresid=resid.maxresid; 
% elapsed time   
  t_elapsed=toc(tinit);
  disp(' ');disp(['Elapsed CPU time is ' num2str(t_elapsed/60) ' minutes'])
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF TwoPhase
%%   internal functions start here. External functions are in 
%%    ./02_functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Molten,phiLAB] = prep_molten(Vark,SurfCoff,Vbgs,zLABid,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extracts the portion of the variable arrays that are below the 
% solidus-liquidus intersection. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        extract phi
         Molten.phi = Vark.phi(zLABid-1:end);         
         Molten.phis = Vark.phis(zLABid-1:end); 
         phiLAB = Vark.phi(zLABid);          
%        adjust the ghosts
         Molten.phi(1)=Molten.phi(2); 
         Molten.phis(1)=Molten.phis(2);          
         Molten.grains = Vark.grains(zLABid-1:end); 
         Molten.grain = Vark.grain(zLABid-1:end);
         Molten.eta=Vark.eta(zLABid-1:end); 
         Molten = permeability(Molten,settings); 
         Molten.SurfTen=SurfCoff(zLABid-1:end); 
         Molten.Vbgs = Vbgs(zLABid-1:end); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = store_molten(Vark,Molten,zLABid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moves everything back to global array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%   zero the global arrays    
    nstag = numel(Vark.phis); 
    n = numel(Vark.phi); 
    Vark.Pex = zeros(n,1); 
    Vark.Ss = zeros(nstag,1); 
    Vark.Vss = zeros(nstag,1);
    Vark.vfs = zeros(nstag,1);
    Vark.Tots = zeros(nstag,1);
            
%   move partially molten section back    
    Vark.Pex(zLABid-1:end)=Molten.Pex; 
    Vark.Ss(zLABid-1:end)=Molten.Ss; 
    Vark.Vss(zLABid-1:end)=Molten.Vss; 
    Vark.vfs(zLABid-1:end)=Molten.vfs; 
    Vark.Tots(zLABid-1:end)=Molten.Tots; 
    
%   calculate the cell-centered velocities
    Vark.S=addghosts(stag(Vark.Ss));Vark.Vs=addghosts(stag(Vark.Vss)); 
    Vark.vf=addghosts(stag(Vark.vfs));Vark.Tot=addghosts(stag(Vark.Ss+Vark.Vss)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dPdz,Sd] = calc_Sd_dPdz(phiL,perms,SurfTen,settings,dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the vertical gradient in excess pressure that enforces the 
% melt flux at the LAB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first, calculate the flux value:
  drhog = settings.g*(settings.rhos - settings.rhof)./settings.mufo*settings.gz;

% now calculate the gradient required to get Sd:   
% generally, S = - k * (dPdz + drhog + dPdz_surfaceTension)
% but in this case, we have things stored so that 
%     S = -k * (dPdz + drhog)  + S_surfacetensioncoefficient * d(phi)/dz
% so just turn that around and solve for dPdz

       
  SsSurf = -0*SurfTen.*(phiL - settings.phimin)/dz; %CJH
  
  Sg = - perms * drhog * (1-phiL).^settings.n1; 
  Sd = Sg * settings.Sd_coefficient ;
 
  dPdz = (Sd - SsSurf - Sg)./(-1*perms)*settings.mufo; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vbgz = upwelling_velocity(Vbgz0,zs,zLAB,DBL_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct vertical solid velocity. This is a trick... only using this to 
% calculate decompression melting, thermal advection terms. Strictly speaking,
% this will violate mass conservation. But maybe that's ok. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   Vbgz = Vbgz0 * 0.5 * (1 + erf((zs - zLAB - 2*DBL_m)./DBL_m)); 
   Vbgz(zs<zLAB) = 0; 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dphidt,dTdt,dCfdt]=SourceTerms(Vark,settings,dCfdt,dphidt,dTdt)
   in_i =2:numel(Vark.Gdot)-1;   
   dphidt.Gdot(in_i)=Vark.Gdot(in_i);
   dTdt.Gdot = zeros(size(Vark.Gdot));
   Melt=1-strcmp(settings.Flags.problem,'One_Phase_T');
   dTdt.Gdot(in_i) = -Melt*settings.L * Vark.Gdot(in_i) ./ Vark.cp(in_i);
   dCfdt.Gdot = zeros(size(Vark.Gdot));
   dCfdt.Gdot(in_i) =(settings.kd_H2O-1)*Vark.Gdot(in_i).*Vark.Cf(in_i)...
                      ./dCfdt.Mix(in_i);                   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dTdt,dphidt,dCfdt]=advective_fluxes(Vark,dz,dt,dCfdt,dphidt,dTdt)
   dTdt.adv = advection_driver(Vark.T,Vark.Tots,dz,dt,'VL_nc');
   dTdt.adv_S = advection_driver(Vark.T,Vark.Ss,dz,dt,'VL_nc');
   dphidt.adv = advection_driver(Vark.phi,Vark.vfs,dz,dt,'VL');
   dCfdt.adv = advection_driver(Vark.Cf,dCfdt.Veff,dz,dt,'VL_nc');
   dCfdt.adv_G = advection_driver(Vark.Cf,dCfdt.Veff-Vark.Vss,dz,dt,'VL_nc');
%    dCfdt.adv_G = advection_driver(Vark.Cf,dCfdt.Veff,dz,dt,'VL_nc');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = MeltingRate(Vark,VTot,dTdt,dCfdt,settings,phimin,dt,z,zLABid,dphidt)
% load in vars of interest
  T = Vark.T; phi = Vark.phi; H2O_f = Vark.Cf; Phy = Vark.P; 
  CO2_f = Vark.Cf_CO2; 
% get solidus and pressure derivatives
  Solidus = SoLiquidus([Phy(2); Phy(2:end)],H2O_f,CO2_f,'katz');
  dT = Solidus.Tlherz - Solidus.Tsol;  
    
%%%%%%%%%%%%%%%%%%%%%%%
% equilibrium melting %
%%%%%%%%%%%%%%%%%%%%%%% 

% calculate denominator 
  dTdH2O=Solidus.dTdH2O;
  if strcmp(settings.Flags.H2O,'variable')
      Hterm = dTdH2O.*H2O_f.*(settings.kd_H2O - 1)./dCfdt.Mix;
  else
      Hterm = 0; 
  end
  Lterm = settings.L ./ (settings.Cp_olv);
  denom = dT + Lterm + Hterm;        
 
% temperature and pressure terms  
  dTdtPsol=-VTot.* settings.rhos .* settings.g .* Solidus.dTdPsol;     
  Tterms = (0*dTdt.adv_S + 0*dTdt.adi+0*dTdt.diff+dTdtPsol); % should be an advective term here?
% water terms
  if strcmp(settings.Flags.H2O,'variable')
      Hterms = -dTdH2O .* (dCfdt.diff+0*dCfdt.adv_G); %CJH
  else
      Hterms = 0; 
  end
% add them up! 
  Gdot = (Tterms + Hterms) ./ denom;  
  
% correct for T < or > Tsol
  Gdot(1:zLABid-1)=0; 
  for iz = zLABid:numel(T)
     if T(iz) < Solidus.Tsol(iz) 
         Gdot(iz) = 0; % crystallization taken care of elsewhere
     elseif T(iz) > Solidus.Tliq(iz)
         if Gdot(iz) < 0 
             Gdot(iz) = 0; 
         end
     elseif T(iz) > Solidus.Tsol(iz)
         if phi(iz) <= settings.phimin && Gdot(iz)<0
             Gdot(iz) = 0;
         end         
     end
  end
  Vark.Gdot = Gdot; 
  Vark.Tsol = Solidus.Tsol; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dt]=dtCFL(vf,dz,cfl,dt_max)
   % calculate time step      
     vfmax = max(abs(vf)); 
     dt = cfl*dz/(vfmax+1e-20); % porous flow max time step
     mindt = dt_max *1e6 * 365 * 24 * 3600; 
     dt = dt*(vfmax~=0)+mindt*(vfmax==0);  % control for vf == 0     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark=permeability(Vark,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           permeability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates cell centered (k) and staggered (ks) permeability using
%    phi     cell centered melt fraction
%    phis    staggered melt fraction
%    grain   cell centered grain size [m]
%    grains  staggered grain size [m]
%    settings.   structure that incluces:
%       np      porosity exponent
%       C       tortuosity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     phi=Vark.phi; phis=Vark.phis; grain=Vark.grain; grains=Vark.grains;
     Vark.perm = (grain.^2).*(phi).^settings.np/settings.C;
     Vark.perms= (grains.^2).*(phis).^settings.np/settings.C;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = darcy(Vark,Vbg,SurfCoff,dz,settings)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               darcy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flow = darcy(phi,phis,ks,P,Vbg,SurfCoff,dz,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns the fluid velocity, darcy flux and solid velocity on cell edges. 
% input
%  phi       melt fraction on cell centers. n+2 elements: first and last
%            elements are ghost cells.
%  phis      melt fraction on cell edges. n+1 elements. phis(1) is at z=0.
%  ks        permeability on cell edges
%  P         excess pressure on cell centers, n+2 elements: first and last
%            elements are ghost cells.
%  Vbg       backround solid velocity on cell edges
%  settings  structure with material properties:
%          .np     permeability exponent. k = phi^np
%          .nxi    compaction viscosity exponent. xi = phi^nxi
%          .n1     switch for 1-phi approximation. (1-phi)^n1
%
% output
%  Flow      structure with fields: 
%      .vfs     fluid velocity on cell edges
%      .Ss      darcy velocity on cell edges
%      .Vss     solid velocity on cell edges
%      .vf      fluid velocity on cell centers
%      .S       darcy velocity on cell centers
%      .Vs      solid velocity on cell centers
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   load in vars of interst
    phi = Vark.phi; phis = Vark.phis; ks = Vark.perms; P = Vark.Pex;     
%   load material parametes    
    n1 = settings.n1; 
    drhog = (settings.rhos - settings.rhof)*settings.g;
    gz = settings.gz;     
%   get useful quantities and derives    
    Fs  = (1-phis).^n1;  
    Fc = (1-phi).^n1;
    P1 =  P.*Fc;
    dPdz = (P1(2:end)-P1(1:end-1))/dz;      
%   calculate velocity of cell faces  
    SsSurf = -SurfCoff.*(phi(2:end)-phi(1:end-1))/dz; 
    %Sd = -settings.Sd_coefficient * phis.^settings.np .* (1-phis).^settings.n1; 
    SPg =-ks.*(dPdz+drhog*Fs*gz); 
    Ss =SPg +SsSurf;            
    Vss = Vbg-Ss; 
    vfs = Ss./(phis+1e-8)+Vss;
%  put it all in the flow structure   
%   staggered values on cell faces
    Vark.Ss = Ss; Vark.Vss = Vss; 
    Vark.vfs = vfs; 
    Vark.Tots = Ss + Vss; 
%   cell centered values
    Vark.S = addghosts(stag(Ss)); Vark.Vs = addghosts(stag(Vss)); 
    Vark.vf = addghosts(stag(vfs)); 
    Vark.Tot = addghosts(stag(Ss + Vss));         
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dTdt,dphidt]=initialize_dT_dphi(nz)
   dTdt.adi = zeros(nz,1); dTdt.diff = zeros(nz,1); 
   dTdt.Gdot = zeros(nz,1); dTdt.adv = zeros(nz,1); 
   dphidt.adv = zeros(nz,1); dphidt.Gdot=zeros(nz,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dCfdt = dCfdt_init(Vark,settings,dz,z) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializes the fluxes for water concentration. Calculates the diffusive 
% fluxes and the src/sink from phase changes. Also calculates the velocity
% that will be used to advect the concentration field later. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% useful quantities:


%   kd = settings.kd_H2O; %CJH
%   kds = kd; 
  
% account for pargasite storage capacity using an effective partitioning
  Cbulk = Vark.phi.*Vark.Cf + Vark.Cs.*(1-Vark.phi); 
  [Vark.Strg,kd]=pargasite_storage(z/1e3,Cbulk,settings,0,'load');
  kds=kd;%(kd(1:end-1)+kd(2:end))/2;     
  dCfdt.kd = kd; 
  dCfdt.kds = kd;   
  phis = Vark.phis; % staggered phi
  phi = Vark.phi; 
  Denom = (phi + kd.*(1-phi)); 
  dCfdt.Mix=Denom; 

  
  
 Vadv = Vark.Vbgzs;
%  Vadv = Vark.Vss;   

% CJH PARG
%   Vadv = Vark.Vbgzs(end); % might want to set this to zero above LAB...
%   
  dCfdt.Veff = (phis.*Vark.vfs + kds.* (1-phis).*Vadv)./(phis + kds.* (1-phis));
   
% % source/sink from phase change (THIS IS CALCULATED LATER)
%   dCfdt.Gdot = Vark.Gdot.*Vark.Cf*(settings.kd_H2O - 1)./dCfdt.Mix;   
  
% calculate diffusive fluxes in fluid and solid  
  Df_eff=settings.Cs_Df./Denom; 
  Ds_eff=settings.Cs_Ds*kd./Denom;
  
  A_f = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,1./Df_eff,1,phi); 
  A_s = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,1./Ds_eff,1,(1-phi)); 
  
  dCfdt_diff_f = (A_f(:,:)) * Vark.Cf(:); % diffusive flux in liquid
  dCfdt_diff_s = (A_s(:,:)) * Vark.Cf(:); % diffusive flux in solid
  
  dCfdt.diff=dCfdt_diff_f+dCfdt_diff_s; 

  dCfdt.maxdiff= max(max([Df_eff.*phi Ds_eff.*(1-phi)])); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vars,Info,kk]=Varstruct(Info,z,outn,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Varstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initializes the variable structure (Vars), output counter (kk) and 
% settings and single-valued time-variables (Info). 
% Add new variables for output here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nz = numel(z); 
        Info = populate_Vars(Info,1,outn,'t','zSOL','zSOLid','zMO','zMOid',...
               'phiLAB','zLABeta','zLABetaid','CbulkResid','Extruded');           
        Info.z = z; 
        Info.z_km = z/1e3;
        
        Vars = populate_Vars([],nz,outn,'S','vf','Vbgz','Pex','Gdot','phi',...
               'eta','P','Vs','T','cp','rho','Kc','dg_um','sig_MPa',...
               'comp','Gu_0_GPa','Tsol','Cs','Cf','Cparg','Cbulk','Strg',...
               'PhiFreeze');         
        kk = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  Vars = populate_Vars(Vars,nz,outn,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% takes list of variables and places them in a structure with initial value
% of zero.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iV = 1:nargin-3
       Vars.(varargin{iV})=zeros(nz,outn);
    end           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vars,Info,kk]=savevar(Vars,Vark,Info,tnow_s,Extruded,k,kk,LABinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              savevar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updates the output structure with current time step values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Vark]=stagger_all(Vark,Vars); 
%   loop through the structure fields and store in the big output structure
    Fields = fieldnames(Vars);
    for iFie = 1:numel(Fields);        
      Vars.(Fields{iFie})(:,kk)= Vark.([Fields{iFie}]);% 's']);
    end    
    Info.t(kk) = tnow_s;
    Info.Extruded(kk) = Extruded; 
    Info.zSOL(kk) = LABinfo.zLAB; Info.zSOLid(kk) = LABinfo.zLABid;
    Info.zLABeta(kk) = LABinfo.zLABeta; Info.zLABetaid(kk) = LABinfo.zLABetaid;
    Info.zMO(kk) = LABinfo.zMO; Info.zMOid(kk) = LABinfo.zMOid;
    Info.phiLAB(kk) = LABinfo.phiLAB;
    
    % check water mass conservation (only bother for output)
    Info.CbulkResid(kk)=0;%calc_CbulkResid(Vark,Vars,Info.z,tnow_s);
    
    kk = kk+1;
    tkyrs=tnow_s/3600/24/365/1e3;
    fprintf('\n   Storing data at step %i, %0.3f kyrs',k,tkyrs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vars,Info]=varfinal(Vars,Info,kk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              varfinal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% removes unfilled cells in output variable arrays. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kk = kk -1;    
    Fields = fieldnames(Vars); 
    for iFie = 1:numel(Fields);        
      Vars.(Fields{iFie})=Vars.(Fields{iFie})(:,1:kk);
    end    
    Info.t = Info.t(1:kk); 
    Info.Extruded = Info.Extruded(1:kk); 
    Info.zSOL = Info.zSOL(1:kk);Info.zSOLid = Info.zSOLid(1:kk);
    Info.zLABeta = Info.zLABeta(1:kk); Info.zLABetaid = Info.zLABetaid(1:kk);
    Info.zMO = Info.zMO(1:kk);Info.zMOid = Info.zMOid(1:kk);
    Info.phiLAB = Info.phiLAB(1:kk);        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = stagger_all(Vark,Vars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         stagger_all -- updates staggered values of all variables that
%                        will be output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fields = fieldnames(Vars);
    for iFie = 1:numel(Fields);
        if isfield(Vars,Fields{iFie})
            Vark.([Fields{iFie} 's'])=stag(Vark.(Fields{iFie}));
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Vark1,resid,tnow_s_1,dt,LABInfo] = split_timestep(Vark,tnow_s,LABInfo,...
%                       settings,nz,z,zs,dz,IVals,BCs,cfl,in_i)
%    dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]
%    phimin =  settings.phimin; 
%    phimax =  settings.phimax; 
%    Vbgs = -settings.Vbg / (1e2 * 365 * 24 * 3600); % [m/s] for mass conservation
%    DBL_m = settings.DBL_m; % mechanical boundary layer thickness at LAB
%    VbgFlag = settings.Flags.VbzFlag;
%    XtalFactor = settings.XtalFactor; 
%    if strcmp(settings.Flags.H2O,'variable'); H2OFlag=1; else H2OFlag = 0; end
% 
% % LAB location info -- need to clean this up.    
% %  the newest                        the old (lag zLAB change for phi)    
%    zLAB1 = LABInfo.zLAB1;            zLAB0 = LABInfo.zLAB0;
%    zLABid1 = LABInfo.zLABid1;        zLABid0 = LABInfo.zLABid0;   
%    phiLAB1 = LABInfo.phiLAB1;        phiLAB0 = LABInfo.phiLAB0;   
%    zLABid = LABInfo.zLABid;  
%    
% %% INITIAL FLUXES   
% %  zero everything 
%    [dTdt,dphidt]=initialize_dT_dphi(nz);
%    
% %  non-advective fluxes      
%    dTdt.adi(in_i) = Vark.Tot(in_i).* settings.dTdz_ad; % adiabatic
%    A_T = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,Vark.rho,Vark.cp,Vark.Kc); 
%    dTdt.diff(:) = (A_T(:,:)) * Vark.T(:); % diffusive
%    dTdt.diff(1) = 0; dTdt.diff(end) = 0; 
%    dCfdt = dCfdt_init(Vark,settings,dz,z); 
%       
% %  calculate stable thermal step
%    dt_TDiff = 0.5*dz*dz/max(Vark.Kc./Vark.rho./Vark.cp);   
%    dt_TAdv = dtCFL(Vark.Tots,dz,cfl,dt_max);
%    dt=min([dt_TDiff dt_TAdv]);
%    
% % calculate advective fluxes
%   [dTdt,dphidt,dCfdt]=advective_fluxes(Vark,dz,dt,dCfdt,dphidt,dTdt);
% 
% % calculate melting/crystallization rate
%    Vgdot = Vark.Vbgz; 
%    Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLAB1,dphidt);   
%    [dphidt,dTdt,dCfdt]=SourceTerms(Vark,settings,dCfdt,dphidt,dTdt);
%    
% % calculate full operators
%    dTdt.Total =  (dTdt.diff + dTdt.adi + dTdt.adv + dTdt.Gdot);
%    plot(Vark.Tsol); hold all; plot(Vark.T)
% %%  update porous flow in sub-steps      
%    if strcmp(settings.Flags.problem,'One_Phase_T')==0
%      tpf = 0; phi0=phimin; 
%      while tpf < dt
% %     zero everything       
%       [dTdt1,dphidt]=initialize_dT_dphi(nz); 
%       dTdt1=dTdt; 
%       dCfdt1 = dCfdt_init(Vark,settings,dz,z); 
%       
% %     calculate stable p.f. step   
%       dt_pf=dtCFL(Vark.vfs,dz,cfl,dt_max); 
%       dt_CfDiff = 0.5*dz*dz/dCfdt.maxdiff;
%       dt_CfCFL=dtCFL(dCfdt.Veff,dz,cfl,dt_max);
%       dt_pf=min([dt_pf dt_CfDiff dt_CfCFL]);
%       if tpf + dt_pf > dt; dt_pf = min([(dt-tpf) dt_pf]); end
%       
% %     advective flux      
%       [dTdt1,dphidt,dCfdt1]=advective_fluxes(Vark,dz,dt_pf,dCfdt1,dphidt,dTdt1);      
% %     source terms
%       Vark = MeltingRate(Vark,Vgdot,dTdt1,dCfdt1,settings,phimin,dt,z,zLAB1,dphidt);
%       [dphidt,dTdt1,dCfdt1]=SourceTerms(Vark,settings,dCfdt1,dphidt,dTdt1);
% %     limit how much phi can change in a step
%       fstep = 0.02;
%       phiproject=abs((dphidt.adv(in_i)*dt + Vark.Gdot(in_i)*dt));
%       if max(abs(phiproject)) > fstep
%           dt_m = min(abs((fstep)./(dphidt.adv(in_i) + Vark.Gdot(in_i))));
%           dt_pf = min([dt_m dt_pf]);
%           if tpf + dt_pf > dt; dt_pf = min([(dt-tpf) dt_pf]); end          
%           [dTdt1,dphidt,dCfdt1]=advective_fluxes(Vark,dz,dt_pf,dCfdt1,dphidt,dTdt1);
%           Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLAB1,dphidt);
%           [dphidt,dTdt1,dCfdt1]=SourceTerms(Vark,settings,dCfdt1,dphidt,dTdt1);
%       end           
% 
% %     update phi and Cf      
%       dphidt.Total = dphidt.adv  + dphidt.Gdot; 
%       dCfdt1.Total = H2OFlag * (dCfdt1.diff + dCfdt1.adv + dCfdt1.Gdot);
%       
%       Old.phi = Vark.phi; 
% 
%       Vark.phi(in_i) = Vark.phi(in_i) + dphidt.Total(in_i)*dt_pf;
%       
%       Old.Cf = Vark.Cf;      
%       
%       Xtal_Gdot=-(Vark.phi-phi0)/dt .* (Vark.T < Vark.Tsol)* XtalFactor*0;      
%       Vark = step_Cf(Vark,Xtal_Gdot,in_i,settings,BCs,dCfdt1,dt_pf*H2OFlag,dz);       
%       Vark = step_xtal_Cf(Vark,settings,dCfdt.kd);      
%       PhiProj=Vark.phi.*(Vark.T>=Vark.Tsol) + phimin * (Vark.T < Vark.Tsol); 
% %       phi0=Vark.phi; 
%                
%       tpf = tpf + dt_pf;
% 
% %     limit phi, set BCs      
%       Vark.phi(Vark.phi<phimin)=phimin;     
%       Vark.phi(Vark.phi>phimax)=phimax;       
%       Vark.phi = BC_setghosts(Vark.phi,BCs.val_phi,BCs.type_phi,dz);
%       
% %     calculate residuals      
%       resid.phi = max(abs(Old.phi(in_i) - Vark.phi(in_i))./Old.phi(in_i));
%       resid.Cf = max(abs(Old.Cf(in_i) - Vark.Cf(in_i))./Old.Cf(in_i));      
%       
%       Vark = update_static_melt(Vark,Vbgs,settings,z,zs,dz,LABInfo,...
%                                 DBL_m,VbgFlag,BCs);
%       dTdt1.Total =  (dTdt1.diff + dTdt1.adi + dTdt1.adv + dTdt1.Gdot); 
% 
%       phinf=Vark.phi(zLABid0-1);
%       if phinf >= Vark.phi(zLABid0) && tpf<dt; dt=tpf; end
%       
%       Vark = MeltingRate(Vark,Vgdot,dTdt1,dCfdt1,settings,phimin,dt,z,zLAB1,dphidt);
%       [LABInfoProject] = find_LAB(Vark,z,settings,LABInfo);
%       if (LABInfoProject.zLABid ~= LABInfo.zLABid0 && tpf<dt) ||...
%               (LABInfoProject.zLABid ~= LABInfo.zLABid && tpf<dt)
%           dt=tpf;
%       end
%       
%       if strcmp(settings.Flags.progress_plot,'yes')
%           Vark.fprog = progress_plot(Vark,z,Vark.fprog);
%       end
%       
%      end
%    
% %    freeze whatever was pushed above the solidus     
%      Xtal_Gdot = -(Vark.phi - phimin)/dt .* (Vark.T < Vark.Tsol)* XtalFactor; 
%      Xtal_L = - settings.L * Xtal_Gdot ./ Vark.cp ;   
%      Vark.phi(1:zLABid0-1)=phimin;
%      Vark.Cbulk= Vark.Cf .* Vark.phi + Vark.Cs .* (1-Vark.phi);
%           
%      Cb=Vark.Cf.*PhiProj + Vark.Cs.*(1-PhiProj);
%      Cb0=settings.Cbulk0;
%      Flux_in=-Vark.Vbgzs(end)*Vark.Cs(end)*(tnow_s+tpf); % CJH PARG
%      Cbz=trapz(z,Cb);
%      Cb_diff=Cbz-trapz(z,Cb0);
%      
%      % Residual: total mass missing relative to current total mass
%      [(Cb_diff - Flux_in)./(Cbz) max(Vark.phi)]
%      % RESID > RESID 2 below --- deifnitely gaining mass up here. 
%    else 
%      Xtal_L = 0; Xtal_Gdot = 0; resid.phi=0;
%      Vark.phi(1:zLABid0-1)=phimin;
%      resid.Cf=0;
%    end      
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% %% update temperature %
% %%%%%%%%%%%%%%%%%%%%%%%
%      
% %   save old values
%     Old.T = Vark.T;
%     
% %   update         
%     dTdtTotal=dTdt.Total + Xtal_L; 
%     Vark.T(in_i) = Vark.T(in_i) + dTdtTotal(in_i)*dt; 
%     Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);
%     Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLAB1,dphidt); 
%      % do this even if not calculating phi evolution to get updated solidus
%     tnow_s_1 = tnow_s + dt;
%     
% %  calculate residuals      
%     resid.T = max(abs(Vark.T(in_i) - Old.T(in_i))./(Old.T(in_i)));        
%     resid.maxresid = max([resid.phi resid.T resid.Cf]);       
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %% calculate "static" quantities %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % density, thermal properties: 
%   Vark = call_MaterialProps(Vark,IVals.Rho_0,IVals.Kc_0,IVals.Cp_0,...
%                             settings,z,dz);  
%   Vark.eta = get_VBR_visc(Vark); 
%   
% % current LAB location
%   [LABInfo] = find_LAB(Vark,z,settings,LABInfo);  
% 
% %  update all, fill the in between with phi  
%    if (LABInfo.zLABid<LABInfo.zLABid0) && strcmp(settings.Flags.problem,'One_Phase_T')==0
% %   get the LAB node   
%     idz=LABInfo.zLABid; % new
%     idz0=LABInfo.zLABid0; % old
%     
% %   account for Cf change with phi-fill 
%     rng=idz:idz0-1; % the range to fill    
%     Vark = step_fill_Cf(Vark,dCfdt.kd,rng,idz0);
%     Vark.Cbulk= Vark.Cf .* Vark.phi + Vark.Cs .* (1-Vark.phi);
% 
% %     Xtal_Gdot=zeros(size(Vark.phi));      
% %     Xtal_Gdot(rng)=(Vark.phi(idz0)-Vark.phi(rng));
%     
%     
% %   now fill it    
%     Vark.phi(rng)=Vark.phi(idz0);
%     
%     
% 
% % %   and adjust Cf appropriately -- IS THIS THE PROBLEM????????  CJH PARG
% %     dCfdt1 = dCfdt_init(Vark,settings,dz,z); % (need the new mixing ratio) 
% %     dCfdt1.Total = 0.0 * (dCfdt1.diff);
% %     Cfid=Vark.Cf(rng);
% %     Vark = step_Cf(Vark,Xtal_Gdot,in_i,settings,BCs,dCfdt1,H2OFlag,dz);
% %     Cfidnew=Vark.Cf(rng);
% %     if H2OFlag == 1; 
% %     disp(['adjusted boundary Cf:' num2str(Cfid') ' to ' num2str(Cfidnew')]);
% %     end
% 
% %   update the LAB nodes
%     LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid;
%    elseif (LABInfo.zLABid>LABInfo.zLABid0) && strcmp(settings.Flags.problem,'One_Phase_T')==0
%        if H2OFlag == 0; 
%         LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid;
%        end
%    end
%    
% %  update other solidus-LAB node infos 
%    LABInfo.zLAB1=LABInfo.zLAB; LABInfo.zLABid1=LABInfo.zLABid;   
%    LABInfo.phiLAB0=LABInfo.phiLAB1; LABInfo.phiLAB1=LABInfo.phiLAB;
% 
% % output the final arrays
%   Vark1 = Vark; 
%   
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vark1,resid,tnow_s_1,dt,LABInfo,Extruded] = ...
       split_timestep(Vark,tnow_s,LABInfo,settings,nz,z,zs,dz,IVals,BCs,...
                      cfl,in_i,Extruded)
   dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]
   phimin =  settings.phimin; 
   XtalFactor = settings.XtalFactor; 
   TempUpdate= 1-strcmp(settings.Flags.TempUpdate,'DiffusionOnly');
   if strcmp(settings.Flags.H2O,'variable'); H2OFlag=1; else H2OFlag = 0; end

% LAB location info -- need to clean this up.    
%  the newest                        the old (lag zLAB change for phi)    
   zLAB1 = LABInfo.zLAB1;            zLAB0 = LABInfo.zLAB0;
   zLABid1 = LABInfo.zLABid1;        zLABid0 = LABInfo.zLABid0;   
   phiLAB1 = LABInfo.phiLAB1;        phiLAB0 = LABInfo.phiLAB0;   
   zLABid = LABInfo.zLABid;  
   
%% INITIAL FLUXES   
%  zero everything 
   [dTdt,dphidt]=initialize_dT_dphi(nz);
   
%  non-advective fluxes      
   dTdt.adi(in_i) = Vark.Tot(in_i).* settings.dTdz_ad; % adiabatic
%    A_T = make_A_matrix_f_ze_CpRhoKc(1/dz/dz,Vark.rho,Vark.cp,Vark.Kc);    
%    dTdt.diff(:) = (A_T(:,:)) * Vark.T(:); % diffusive
   dTdt.diff = make_dTdt_diff(z,Vark.rho,Vark.cp,Vark.Kc,Vark.T);% * TempUpdate; % CJH
   dTdt.diff(1) = 0; dTdt.diff(end) = 0; 
   dCfdt = dCfdt_init(Vark,settings,dz,z); 
      
%  calculate stable thermal step
   dt_TDiff = 0.5*dz*dz/max(Vark.Kc./Vark.rho./Vark.cp);   
   dt_TAdv = dtCFL(Vark.Tots,dz,cfl,dt_max);
   dt_adi = 0.1 * min(abs(Vark.T./dTdt.adi));
   dt=min([dt_TDiff dt_TAdv dt_adi]);
   
% calculate advective fluxes
  [dTdt,dphidt,dCfdt]=advective_fluxes(Vark,dz,dt,dCfdt,dphidt,dTdt);
   
% calculate melting/crystallization rate
   Vgdot = Vark.Vbgz; 
   Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLABid0,dphidt);   
   [dphidt,dTdt,dCfdt]=SourceTerms(Vark,settings,dCfdt,dphidt,dTdt);
   
% calculate full operators
   dTdt.Total =  (dTdt.diff + dTdt.adi*TempUpdate + dTdt.adv*TempUpdate + dTdt.Gdot*TempUpdate);
%    plot(Vark.Tsol); hold all; plot(Vark.T)
%%  update porous flow in sub-steps   
   
   if strcmp(settings.Flags.problem,'One_Phase_T')==0 && TempUpdate == 1
       
%    step forward phi and Cf       
     [Vark,resid,tpf,dt,LABInfo,Xtal_L,Extruded] = ...
         timestep_phi_Cf(Vark,settings,LABInfo,BCs,dTdt,dCfdt,cfl,...
                         in_i,dt,dz,nz,z,zs,Extruded);
   
          
%      Cb=Vark.Cbulk;
%      Cb0=settings.Cbulk0;
%      Flux_in=-Vark.Vbgzs(end)*Vark.Cs(end)*(tnow_s+tpf); % CJH PARG
%      Cbz=trapz(z,Cb);
%      Cb_diff=Cbz-trapz(z,Cb0);
     
     % Residual: total mass missing relative to current total mass
%      [(Cb_diff - Flux_in)./(Cbz) max(Vark.phi)]
     % RESID > RESID 2 below --- deifnitely gaining mass up here. 
   else 
     Xtal_L = 0; Xtal_Gdot = 0; resid.phi=0;
     Vark.phi(1:zLABid0-1)=phimin;
     resid.Cf=0;
   end      

%%%%%%%%%%%%%%%%%%%%%%%
%% update temperature %
%%%%%%%%%%%%%%%%%%%%%%%
     
%   save old values
    Old.T = Vark.T;
    
%   update             
    dTdtTotal=(dTdt.Total + Xtal_L*TempUpdate) ; 
    Vark.T(in_i) = Vark.T(in_i) + dTdtTotal(in_i)*dt; 
    Vark.T = BC_setghosts(Vark.T,BCs.val_T,BCs.type_T,dz);
    Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLABid0,dphidt); 
     % do this even if not calculating phi evolution to get updated solidus
    tnow_s_1 = tnow_s + dt;
    
%  calculate residuals      
    resid.T = max(abs(Vark.T(in_i) - Old.T(in_i))./(Old.T(in_i)));        
    resid.maxresid = max([resid.phi resid.T resid.Cf]);       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% calculate "static" quantities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% density, thermal properties: 
  Vark = call_MaterialProps(Vark,IVals.Rho_0,IVals.Kc_0,IVals.Cp_0,...
                            settings,z,dz);  
  Vark.eta = get_VBR_visc(Vark); 
  
% current LAB location
  [LABInfo] = find_LAB(Vark,z,settings,LABInfo);  
  
%  update all, fill the in between with phi  
   if (LABInfo.zLABid0-LABInfo.zLABid>1) && strcmp(settings.Flags.problem,'One_Phase_T')==0
      LABInfo.zLABid=LABInfo.zLABid0-1;
      LABInfo.zLAB=z(LABInfo.zLABid);
   end
   if (LABInfo.zLABid<LABInfo.zLABid0) && strcmp(settings.Flags.problem,'One_Phase_T')==0
%   get the LAB node   
    idz=LABInfo.zLABid; % new
    idz0=LABInfo.zLABid0; % old
   
%     Xtal_Gdot=zeros(size(Vark.phi));      
%     Xtal_Gdot(rng)=(Vark.phi(idz0)-Vark.phi(rng));

%   only fill if latent heat is accounted for     
    phifill = Vark.phi(idz0); % the phi to fill with 
    TLab = Vark.T(idz);     
    dTlatent = phifill * settings.L / Vark.cp(idz); % ???
    Tsolz = Vark.Tsol(idz); 
    Tcrit = (TLab - dTlatent) - Tsolz; 
%     disp(Tcrit)
    Tcrit = 10; 
    if Tcrit > 0   % then move it 
%         disp('Tcrit > 0, move LAB')        
%      now fill it    
       Vark.phi(idz)=Vark.phi(idz0);     
       
       % %   and adjust Cf appropriately -- IS THIS THE PROBLEM????????  CJH PARG
       %     dCfdt1 = dCfdt_init(Vark,settings,dz,z); % (need the new mixing ratio)
       %     dCfdt1.Total = 0.0 * (dCfdt1.diff);
       %     Cfid=Vark.Cf(rng);
       %     Vark = step_Cf(Vark,Xtal_Gdot,in_i,settings,BCs,dCfdt1,H2OFlag,dz);
       %     Cfidnew=Vark.Cf(rng);
       %     if H2OFlag == 1;
       %     disp(['adjusted boundary Cf:' num2str(Cfid') ' to ' num2str(Cfidnew')]);
       %     end

%       update the LAB nodes
        LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid;
    else
        LABInfo.zLAB = LABInfo.zLAB0; 
        LABInfo.zLABid = LABInfo.zLABid0; 
    end
    
   elseif (LABInfo.zLABid>LABInfo.zLABid0) && strcmp(settings.Flags.problem,'One_Phase_T')==0
       if H2OFlag == 0; 
        LABInfo.zLAB0=LABInfo.zLAB; LABInfo.zLABid0=LABInfo.zLABid;
       end
   end
   
%  update other solidus-LAB node infos 
   LABInfo.zLAB1=LABInfo.zLAB; LABInfo.zLABid1=LABInfo.zLABid;   
   LABInfo.phiLAB0=LABInfo.phiLAB1; LABInfo.phiLAB1=LABInfo.phiLAB;

% output the final arrays
  Vark1 = Vark; 
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  Vark = update_static_melt(Vark,Vbgs,settings,z,zs,dz,LABInfo,...
                                DBL_m,VbgFlag,BCs)                                
%    LAB Info				
     zLAB=LABInfo.zLAB0;zLABid=LABInfo.zLABid0; phiLAB=LABInfo.phiLAB0;
     if strcmp(settings.Flags.LABdef,'visc')==1
         zLAB_Vbg=LABInfo.zLABeta;
     elseif strcmp(settings.Flags.LABdef,'phi')==1
         zLAB_Vbg=LABInfo.zLAB;
     end
% melt system:   
%    calculate pemeability        
     Vark.phis = stag_phi(Vark,settings);
     Vark=permeability(Vark,settings);
%    calculate surface tension coefficient on staggered grid
     SurfCoff=SurfaceTension(Vark,settings); 

     if strcmp(VbgFlag,'var_z_con_t')==0
         if strcmp(VbgFlag,'variable')
             Vark.Vbgzs=upwelling_velocity(Vbgs,zs,zLAB_Vbg,DBL_m);
             Vark.Vbgz = addghosts(stag(Vark.Vbgzs));
         elseif strcmp(VbgFlag,'con_asthen')
             Vark.Vbgzs = Vbgs.*ones(size(Vark.phis));
             Vark.Vbgzs(zs<zLAB_Vbg)=0;
             Vark.Vbgz = addghosts(stag(Vark.Vbgzs));
         end
     end
     
     if strcmp(settings.Flags.problem,'Two_Phase_Sd')
%       extract only the section below the solidus-liquidus intersection
        [Molten,phiLAB] = prep_molten(Vark,SurfCoff,Vark.Vbgzs,zLABid,settings);
%       calculate the boundary conditions!          
        Sd_type = BCs.type_P; Sd_val = BCs.val_P; 
        Sd_type(1) = 2; % it's a flux!
        [Sd_val(1),Sdflux] = calc_Sd_dPdz(Molten.phi(2),Molten.perms(1),...
                    Molten.SurfTen(2),settings,dz);          
%       and now do the calculation         
        [Molten.Pex,A,src] = P_solve_gen(Molten.phi,Molten.phis,Molten.perms,...
                     Molten.eta,Molten.SurfTen,dz,settings,Sd_val,Sd_type);
%       Get flux and velocities    
        Molten = darcy(Molten,Molten.Vbgs,Molten.SurfTen,dz,settings);          
%       store in global array
        Vark = store_molten(Vark,Molten,zLABid);         
     else
        [Vark.Pex,A,src] = P_solve_gen(Vark.phi,Vark.phis,Vark.perms,...
                               Vark.eta,SurfCoff,dz,settings,BCs.val_P,BCs.type_P); 
%       Get flux and velocities    
        Vark = darcy(Vark,Vark.Vbgzs,SurfCoff,dz,settings);      
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = step_Cf(Vark,Xtal_Gdot,in_i,settings,BCs,dCfdt,dt,dz)

%   Cbulk0= Vark.Cf .* Vark.phi + Vark.Cs .* (1-Vark.phi);
  kd=dCfdt.kd;
%   Xtal_Cf =(kd-1).*Xtal_Gdot.*Vark.Cf./dCfdt.Mix;
%   Vark.Cf(in_i) = Vark.Cf(in_i) + (dCfdt.Total(in_i)+Xtal_Cf(in_i))*dt;  
  Vark.Cf(in_i) = Vark.Cf(in_i) + dCfdt.Total(in_i)*dt;  
  Vark.Cs = Vark.Cf.*kd; 
%   Cbulk1= Vark.Cf .* Vark.phi + Vark.Cs .* (1-Vark.phi);
  
% fill up amphibole storage -- all water mass goes into pargasite
% until pargasite storage capacity is exceeded. 
  Cbulk= Vark.Cf .* Vark.phi + Vark.Cs .* (1-Vark.phi);
  Vark.Cparg = Vark.Cparg + Cbulk.* (Vark.Strg>Vark.Cparg);
  Vark.Cf(Vark.Strg>Vark.Cparg)=0; 

  
  Vark.Cs = Vark.Cf.*kd;
  Vark.Cs = BC_setghosts(Vark.Cs,BCs.val_Cs,BCs.type_Cs,dz);
  Vark.Cf=BC_setghosts(Vark.Cf,BCs.val_Cs./settings.kd_H2O,BCs.type_Cs,dz);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = step_xtal_Cf(Vark,phi1,phi2,settings,kdz)
% for crystallization when T < Tsol

% % get phi info
%   phi1 = Vark.phi; 
%   phimin = settings.phimin; 

% loop through, conserve bulk water when crystallizing  
  Cf0=Vark.Cf; 
  Cf1=Vark.Cf;
  for iz = 1:numel(phi1)     
     
     if Vark.T(iz) < Vark.Tsol(iz)
%           get partition coefficient         
            if numel(kdz) > 1; kd = kdz(iz); else kd = kdz; end
%           conserve bulk water in the crystallization step         
            Bulk = phi1(iz).*Vark.Cf(iz) + (1-phi1(iz)).*Vark.Cf(iz)*kd; 
%           crystallize to phimin
            
            Vark.Cf(iz) = Bulk ./ (phi2(iz) + (1-phi2(iz))*kd);         
            Cf1(iz)=Vark.Cf(iz);
            Vark.Cs(iz) = Vark.Cf(iz) * kd;
     end
  end
  


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = step_fill_Cf(Vark,kdz,rng,idz)  
% for adjusting Cf when LAB shallows

% get phi info
  phi1 = Vark.phi; 
  phifill=Vark.phi(idz); % the value to put in 

% loop through, conserve bulk water when crystallizing  
  for iz = rng(1):rng(end)                
%     get partition coefficient         
      if numel(kdz) > 1; kd = kdz(iz); else kd = kdz; end
%     conserve bulk water: initial value
      Bulk = phi1(iz).*Vark.Cf(iz) + (1-phi1(iz)).*Vark.Cf(iz)*kd; 
%     fill to phifill
      Vark.Cf(iz) = Bulk ./ (phifill + (1-phifill)*kd); 
      Vark.Cs(iz) = Vark.Cf(iz) * kd;     
  end
  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vark,resid,tpf,dt,LABInfo,Xtal_L,Extruded] = ...
    timestep_phi_Cf(Vark,settings,LABInfo,BCs,dTdt,dCfdt,cfl,...
                    in_i,dt,dz,nz,z,zs,Extruded)

% LAB location info 
zLAB1 = LABInfo.zLAB1;
zLABid0 = LABInfo.zLABid0;

% flags and stuff 
  phimin=settings.phimin;phimax=settings.phimax;
  dt_max = settings.dt_max; % max step for advection if velocity = 0 [Myrs]
  Vgdot = Vark.Vbgz;
  XtalFactor = settings.XtalFactor;
  Vbgs = -settings.Vbg / (1e2 * 365 * 24 * 3600); % [m/s] for mass conservation
  DBL_m = settings.DBL_m; % mechanical boundary layer thickness at LAB
  VbgFlag = settings.Flags.VbzFlag;

% the time stepping  
tpf = 0;  keepgoing='gogogo';
while strcmp(keepgoing,'gogogo') || strcmp(keepgoing,'recalc')
%   save initial values
    Init.Vark = Vark;
    Init.LABInfo=LABInfo; 
    Init.tpf= tpf; 
    
%   initialize the derivative structures
    [dTdt1,dphidt]=initialize_dT_dphi(nz);
    dTdt1=dTdt;
    dCfdt1 = dCfdt_init(Vark,settings,dz,z);
    
%   calculate stable p.f. step
    dt_pf=dtCFL(Vark.vfs,dz,cfl,dt_max);
    if tpf + dt_pf > dt; dt_pf = min([(dt-tpf) dt_pf]); end

%   advective flux
    [dTdt1,dphidt,dCfdt1]=advective_fluxes(Vark,dz,dt_pf,dCfdt1,dphidt,dTdt1);
    
%   source terms
    Vark = MeltingRate(Vark,Vgdot,dTdt1,dCfdt1,settings,phimin,dt,z,zLABid0,dphidt);
    [dphidt,dTdt1,dCfdt1]=SourceTerms(Vark,settings,dCfdt1,dphidt,dTdt1);
    
%   limit how much phi can change in a step
    fstep = 0.02;
    phiproject=abs((dphidt.adv(in_i)*dt + Vark.Gdot(in_i)*dt));
    if max(abs(phiproject)) > fstep
        dt_m = min(abs((fstep)./(dphidt.adv(in_i) + Vark.Gdot(in_i))));
        dt_pf = min([dt_m dt_pf]);
        if tpf + dt_pf > dt; dt_pf = min([(dt-tpf) dt_pf]); end
        [dTdt1,dphidt,dCfdt1]=advective_fluxes(Vark,dz,dt_pf,dCfdt1,dphidt,dTdt1);
        Vark = MeltingRate(Vark,Vgdot,dTdt,dCfdt,settings,phimin,dt,z,zLABid0,dphidt);
        [dphidt,dTdt1,dCfdt1]=SourceTerms(Vark,settings,dCfdt1,dphidt,dTdt1);
    end
    
%   update phi 
    dphidt.Total = dphidt.adv  + dphidt.Gdot;
    PartMolt=zLABid0-1:nz;
    Vark.phi(PartMolt) = Vark.phi(PartMolt) + dphidt.Total(PartMolt)*dt_pf;
%     Vark.phi(in_i) = Vark.phi(in_i) + dphidt.Total(in_i)*dt_pf;
    
%   update time    
    tpf = tpf + dt_pf;
    
%   limit phi, set BCs
    for iz = 1:nz
        phiz=Vark.phi(iz);         
        if Vark.T(iz) > Vark.Tsol(iz)
           Vark.phi(iz) = phiz * (phiz > 1e-4) + 1e-4 * (phiz <= 1e-4);  
        end        
    end
    Vark.phi(Vark.phi<phimin)=phimin;
    Vark.phi(Vark.phi>phimax)=phimax;
    
    Vark.phi = BC_setghosts(Vark.phi,BCs.val_phi,BCs.type_phi,dz);
    
%   calculate residuals
    resid.phi = max(abs(Init.Vark.phi(in_i) - Vark.phi(in_i))./Init.Vark.phi(in_i));
    resid.Cf = 0;
    
%   update two-phase variables     
    Vark = update_static_melt(Vark,Vbgs,settings,z,zs,dz,LABInfo,...
        DBL_m,VbgFlag,BCs);
            
%   check how large infilitration phi is    
    phinf=Vark.phi(zLABid0-1);
    if phinf >= Vark.phi(zLABid0) && tpf<dt; 
        dt=tpf; 
        keepgoing='stop';
    end
                 
%   check if max time has elapsed
    if tpf >= dt;
        keepgoing='stop';
    end
    
%   plot progress (if desired)    
    if strcmp(settings.Flags.progress_plot,'yes')
        Vark.fprog = progress_plot(Vark,z,Vark.fprog);
    end   
    
end    

% freeze whatever was pushed above the solidus
%   calculate crystallization, heating rates  
    XtalF=calc_XtalFactor(z,settings);
    Xtal_Gdot = zeros(size(Vark.T)); 
    Xtal_Gdot(2:zLABid0-1) = -(Vark.phi(2:zLABid0-1)-phimin)/dt ...
                             .* XtalF(2:zLABid0-1);
    Vark.PhiFreeze(2:zLABid0-1) = Vark.PhiFreeze(2:zLABid0-1) + (Vark.phi(2:zLABid0-1)-phimin)...
                             .* XtalF(2:zLABid0-1);
    Extruded = Extruded + sum((1 - XtalF(2:zLABid0-1)).*Vark.phi(2:zLABid0-1)*dz);
               % units of length [m]. i.e., thickness of intruded sheet. 
               % total volume would be (areal extent of eruption * thickness). 
    % intrusion length scale - how many nodes to distibute heat over?
      intru_length=settings.DikingBL_km;
      dz = (z(2)-z(1))/1000; 
      Nz = round(intru_length/dz); 
      Total_Energy = sum(- settings.L * Xtal_Gdot ./ Vark.cp  * dz); 
      Heat_src = Total_Energy / intru_length; 
      Xtal_L = zeros(size(z)); 
      if zLABid0-Nz-1 > 1
          startid= zLABid0-Nz-1; 
      else
          startid = 2; 
      end
      Xtal_L(startid:zLABid0-1) = Heat_src;
    
%     Xtal_L = - settings.L * Xtal_Gdot ./ Vark.cp ;
    Vark.phi(1:zLABid0-1)=phimin;    
%     disp(num2str([max(Xtal_L) phinf]))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Resid=calc_CbulkResid(Vark,Vars,z,t)
  Cb=Vark.Cbulk;
  Cb0=Vars.Cbulk(:,1); 
  Cparg=Vark.Cparg; 
  Cparg0=Vars.Cparg(:,1); 
  
%   Flux_in=-Vark.Vs(end)*Vark.Cs(end)*t; % end values don't change through time
  Flux_in=-Vark.Vbgzs(end)*Vark.Cs(end)*t; % CJH PARG
  
  P1=trapz(z,Cparg);
  P2=trapz(z,Cparg0); 
  Parg_Sink=P1-P2;
  
  
  Cbz=trapz(z,Cb);
  Cb_diff=Cbz-trapz(z,Cb0);

  % Residual: total mass missing relative to current total mass
  Resid = (Cb_diff - Flux_in + Parg_Sink)./(Cbz+P1); 
end

function [keepgoing,Info] = check_stop_run(keepgoing,Info,resid,ss_tol,tnow_s,...
                                    t_max_Myrs,LABInfo,settings,Vark,kt_test)
    
    t_s_to_Myrs = 3600*24*365*1e6; % seconds --> Myr factor
    if resid.maxresid<=ss_tol; 
        Info.final_message='Reached Steady State';                 
        keepgoing = 0; 
    end
    if tnow_s/t_s_to_Myrs > t_max_Myrs
       keepgoing = 0;        
       Info.final_message=['Reached desired t: ' num2str(tnow_s/t_s_to_Myrs) '[Myr]'];
    end    
    if LABInfo.zLABeta/1e3 <= settings.Z_moho_km || ...
        LABInfo.zLAB/1e3 <= settings.Z_moho_km
        keepgoing = 0;         
        Info.final_message=['thinned to Moho, quitting']; 
    elseif LABInfo.zLABeta/1e3 <= settings.z_thin || ...
            LABInfo.zLAB/1e3 <= settings.z_thin
        keepgoing = 0;        
        Info.final_message=['thinned to ' num2str(settings.z_thin) ' km, quitting']; 
    end
    if max(Vark.phi)>=settings.phimax && strcmp(settings.Flags.phikill,'stop_on_max_phi')
        keepgoing = 0;         
        Info.final_message='hit max phi';
    end
    if LABInfo.zLABeta/1e3 >= settings.Z_LAB_max
       keepgoing=0;        
       Info.final_message='zLAB deepened to target';
    end
    if strcmp(settings.Flags.problem(1:9),'Two_Phase')
        if sum(isnan(Vark.phi))>0                        
            Info.final_message='NAN ATTACK'; 
            keepgoing = 0;
        end
    end
    if kt_test == 0 
        Info.final_message='reached maximum time steps';
    end
    if strcmp(settings.Flags.problem,'One_Phase_T')==0
        Molten=sum(Vark.T>Vark.Tsol);
        if Molten<=3
            disp('T < Tsol everywhere! Cannot calculate melt migration')
            disp('switching to single phase T evolution...')
            keepgoing = 1; % exit for time step loop  
            settings.Flags.problem='One_Phase_T';
        end                
    end
    if strcmp(settings.Flags.problem,'One_Phase_T')==1 && ...
       strcmp(settings.Flags.problemkill,'stop_on_no_melt')==1
            keepgoing = 0; % exit for time step loop              
            Info.final_message='T < Tsol everywhere, quitting';
    end
end
