%% ===================================================================== %%
%%                     CB_001_0D_scalar.m
%% ===================================================================== %%
%  Calls VBR using a single thermodynamic state
%% ===================================================================== %%
   clear

%% ====================================================
%% Load and set VBR parameters ========================
%% ====================================================

% put VBR in the path
  path_to_top_level_vbr='../../';
  addpath(path_to_top_level_vbr)
  vbr_init

%  write method list (these are the things to calculate)
%  all methods will end up as output like:
%      VBR.out.elastic.anharmonic, VBR.out.anelastic.eBurgers, etc.
   VBR.in.elastic.methods_list={'anharmonic';'poro_Takei';'SLB2005'};
   VBR.in.viscous.methods_list={'HK2003'};
   VBR.in.anelastic.methods_list={'AndradePsP';'YT_maxwell'};

%  load anharmonic parameters, adjust Gu_0_ol
%  all paramss in ../4_VBR/VBR_version/params/ will be loaded in call to VBR spine,
%  but you can load them here and adjust any one of them (rather than changing those
%  parameter files).
   VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
   VBR.in.elastic.anharmonic.Gu_0_ol = 75.5; % olivine reference shear modulus [GPa]

%  frequencies to calculate at
   VBR.in.SV.f = logspace(-2.2,-1.3,4);

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

% HF settings
  HF.Tsurf_C=0; % surface temperature [C]
  HF.Tasth_C=1350; % asthenosphere temperature [C]
  HF.V_cmyr=8; % half spreading rate [cm/yr]
  HF.Kappa=1e-6; % thermal diffusivity [m^2/s]
  HF.rho=3300; % density [kg/m3]
  HF.t_Myr=linspace(0,500,100)+1e-12; % seaflor age [Myrs]
  HF.z_km=linspace(0,200,50)'; % depth, opposite vector orientation [km]

% HF calculations
  HF.s_in_yr=(3600*24*365); % seconds in a year [s]
  HF.t_s=HF.t_Myr*1e6*HF.s_in_yr; % plate age [s]
  HF.x_km=HF.t_s / (HF.V_cmyr / HF.s_in_yr / 100) / 1000; % distance from ridge [km]

  % calculate HF cooling model for each plate age
  HF.dT=HF.Tasth_C-HF.Tsurf_C;
  HF.T_C=zeros(numel(HF.z_km),numel(HF.x_km));
  for HF.i_t = 1:numel(HF.t_s)
    HF.erf_arg=HF.z_km*1000/(2*sqrt(HF.Kappa*HF.t_s(HF.i_t)));
    HF.T_C(:,HF.i_t)=HF.Tsurf_C+HF.dT * erf(HF.erf_arg);
  end


% store in VBR state variables
  % set HF temperature, convert to K
  VBR.in.SV.T_K = HF.T_C+273;
  % construct pressure
  HF.P_z=HF.rho*9.8*HF.z_km*1e3/1e9; %
  VBR.in.SV.P_GPa = repmat(HF.P_z,1,numel(HF.t_s)); % pressure [GPa]

% set the other state variables (ISV) as scalar matrices of same size
  sz=size(HF.T_C);
  VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
  VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
  VBR.in.SV.chi=1*ones(sz); % composition fraction: 1 for olivine, 0 for crust
  VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction
  VBR.in.SV.dg_um = 0.01 * 1e6 * ones(sz); % grain size [um]
  VBR.in.SV.Ch2o = 0 * ones(sz) ; % water concentration

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   [VBR] = VBR_spine(VBR) ;

%% ====================================================
%% Display some things ================================
%% ====================================================
  figure()
  ax1=subplot(2,2,1)
  contourf(HF.t_Myr,HF.z_km,HF.T_C,20)
  colormap(ax1,hot)
  xlabel('Seaflor Age [Myr]')
  ylabel('Depth [km]')
  set(gca,'ydir','rev')
  title('Temperature [C]')
  colorbar()

  for i_f=1:3
     ax=subplot(2,2,i_f+1)
     contourf(HF.t_Myr,HF.z_km,VBR.out.anelastic.AndradePsP.V(:,:,i_f)/1e3,20)
     colormap(ax,winter);
     xlabel('Seaflor Age [Myr]')
     ylabel('Depth [km]')
     set(gca,'ydir','rev')
     title(['V_s [km/s] AndradePsP at ',num2str(VBR.in.SV.f(i_f)),' Hz'])
     colorbar()
  end

  dV=abs(VBR.out.anelastic.AndradePsP.V-VBR.out.anelastic.YT_maxwell.V);
  dV=dV./VBR.out.anelastic.YT_maxwell.V*100;
  figure()
  for i_f=1:4
     subplot(2,2,i_f)
     contourf(HF.t_Myr,HF.z_km,log10(dV(:,:,i_f)+1e-22),50)
     colormap(winter)
     xlabel('Seaflor Age [Myr]')
     ylabel('Depth [km]')
     set(gca,'ydir','rev')
     title(['log_1_0(perc. diff.) between Andrade, Maxwell at ',num2str(VBR.in.SV.f(i_f)),' Hz'])
     colorbar()
  end
