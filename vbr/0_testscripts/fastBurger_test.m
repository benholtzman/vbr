%% ===================================================================== %%
%%                     CB_003_JF10.m
%% ===================================================================== %%
%  Calls VBR after setting the elastic modulus to match the JF10
%  modulus.
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
   VBR.in.elastic.methods_list={'anharmonic'};
   VBR.in.anelastic.methods_list={'eBurgers'};
   VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
   VBR.in.anelastic.eBurgers.method='FastBurger';

%  frequencies to calculate at
   VBR.in.SV.f = 1./logspace(-2,4,100);

%% ====================================================
%% Define the Thermodynamic State =====================
%% ====================================================

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   VBR.in.SV.T_K=700:50:1200;
   VBR.in.SV.T_K=VBR.in.SV.T_K+273;
   sz=size(VBR.in.SV.T_K); % temperature [K]

%  intensive state variables (ISV)
   VBR.in.SV.dg_um=3.1*ones(sz);
   VBR.in.SV.P_GPa = 0.2 * ones(sz); % pressure [GPa]
   VBR.in.SV.rho = 3300 * ones(sz); % density [kg m^-3]
   VBR.in.SV.sig_MPa = 10 * ones(sz); % differential stress [MPa]
   VBR.in.SV.chi=ones(sz); % composition fraction  1 = olivine, 0 = crust

%  structural state variables (SSV)
   VBR.in.SV.phi = 0.0 * ones(sz); % melt fraction

%  compositional state variables (CSV)
   VBR.in.SV.Ch2o = 0 * ones(sz) ; % water concentration

%% ====================================================
%% CALL THE VBR CALCULATOR ============================
%% ====================================================

   % run it initially (eBurgers uses high-temp background only by default)
   [VBR_fastBurger] = VBR_spine(VBR) ;

   % % adjust VBR input and get out eBurgers with background + peak
   VBR.in.anelastic.eBurgers.method='PointWise';
   [VBR_slowBurger] = VBR_spine(VBR) ;

%% ====================================================
%% Display some things ================================
%% ====================================================

  diff1=abs(VBR_fastBurger.out.anelastic.eBurgers.Qinv-VBR_slowBurger.out.anelastic.eBurgers.Qinv);
  disp('max Qinv difference, abs(fast - slow):')
  disp(max(diff1(1,:)))
  disp('max Qinv % difference, 100* abs(fast - slow) / slow:')
  smallQ=1e-4;
  diff1=diff1(:)./(VBR_slowBurger.out.anelastic.eBurgers.Qinv+smallQ)*100;
  disp(max(diff1(1,:)))

  figure
  subplot(1,2,1)
  loglog(VBR.in.SV.f,squeeze(VBR_fastBurger.out.anelastic.eBurgers.Qinv(1,6,:)),'k','displayname','Fast')
  hold all
  loglog(VBR.in.SV.f,squeeze(VBR_slowBurger.out.anelastic.eBurgers.Qinv(1,6,:)),'--r','displayname','Slow')
xlabel('f [Hz]'); ylabel('Qinv [GPa]')

  subplot(1,2,2)
  loglog(VBR.in.SV.f,squeeze(VBR_fastBurger.out.anelastic.eBurgers.M(1,6,:)),'k','displayname','Fast')
  hold all
  loglog(VBR.in.SV.f,squeeze(VBR_slowBurger.out.anelastic.eBurgers.M(1,6,:)),'--r','displayname','Slow')
  legend('location','SouthEast')
  xlabel('f [Hz]'); ylabel('M [GPa]')
