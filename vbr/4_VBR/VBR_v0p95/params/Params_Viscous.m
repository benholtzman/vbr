%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Params_visc.m
%
% input is optional with phi_c and x_phi_c as first and second optional 
% inputs. 
%
% x_phi_c: correction for nominally melt free lab experiements. A value
%           of speefac should be >= 1. A value of 1 does not change the
%           flow law, a value > 1 will decrease strain rate, accounting
%           for experimental conditions that were not trully melt free
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = Params_Viscous(method)

% hirth and kohlstedt 2003
if strcmp(method,'HK2003') 
%  set melt effects    
   phi_c = [1e-5 1e-5 1e-5]; % diff disl gbs
   x_phi_c = [5 1 5/2]; 
   
%  load standard constants
   params = load_HK03_flowlaw_constants(phi_c,x_phi_c); 
   
%  other settings   
   params.ch2o_o = 50; % reference water content [ppm] ("dry" below this value)
   params.P_dep_calc='yes'; % pressure-dependent calculation? 'yes' or 'no'.      
end

% hansen et al., 
if strcmp(method,'LH2012')
   %  set melt effects    
   phi_c = [0.001 0.001 0.001];
   x_phi_c = [5 1 5/2]; 
   
%  load standard constants
   params = load_LH12_flowlaw_constants(phi_c,x_phi_c); 
   
%  other settings      
   params.P_dep_calc='yes'; % pressure-dependent calculation? 'yes' or 'no'.      
end


end


function params = load_HK03_flowlaw_constants(phi_c,x_phi_c)

%% Coble Diffusion creep (GB)
%  dry
   params.diff.A = 1.5e9 ; % preexponential for coble diffusion creep
   params.diff.Q = 375e3 ;% activation energy for coble diffusion creep
   params.diff.V = 10e-6 ; % activation volume for coble diff creep 
   params.diff.p = 3 ; % grain size exponent
   params.diff.n = 1 ; % stress exponent
   params.diff.r = 0 ; % water fugacity exponent
   params.diff.alf = 25 ; % melt factor
   params.diff.phi_c=phi_c(1);
   params.diff.x_phi_c=x_phi_c(1);
   
%  wet   
   params.diff.A_wet = 2.5e7 ; % preexponential for coble diffusion creep
   params.diff.Q_wet = 375e3 ;% activation energy for coble diffusion creep
   params.diff.V_wet = 10e-6 ; % activation volume for coble diff creep
   params.diff.p_wet = 3 ; % grain size exponent
   params.diff.n_wet = 1 ; % stress exponent   
   params.diff.r_wet = 0.7 ; % water fugacity exponent
   params.diff.alf_wet = 25 ; % melt factor
   params.diff.phi_c_wet=phi_c(1);
   params.diff.x_phi_c_wet=x_phi_c(1);
   
%% Dislocation creep 
%  dry
   params.disl.A = 1.1e5 ; % preexponential  
   params.disl.Q = 530e3 ;% activation energy
   params.disl.V = 15e-6 ; % activation volume (HK03 doesn't report, using value from LH12 here)
   params.disl.n = 3.5 ; % stress exponent
   params.disl.p = 0; % grain size exponent
   params.disl.alf = 30 ; % melt factor
   params.disl.r = 0 ; % water fugacity exponent
   params.disl.phi_c=phi_c(2);
   params.disl.x_phi_c=x_phi_c(2);
%  wet
   params.disl.A_wet = 1600.0 ; % preexponential
   params.disl.Q_wet = 520e3 ;% activation energy
   params.disl.V_wet = 22e-6 ; % activation volume
   params.disl.n_wet = 3.5 ; % stress exponent
   params.disl.p_wet = 0; % grain size exponent
   params.disl.alf_wet = 30 ; % melt factor
   params.disl.r_wet = 1.2 ; % water fugacity exponent
   params.disl.phi_c_wet=phi_c(2);
   params.disl.x_phi_c_wet=x_phi_c(2);
   
%% GBS disl accomodated (dry only)
   params.gbs.A_lt1250 = 6500  ; % preexponential for GBS-disl creep
   params.gbs.Q_lt1250 = 400e3  ; % activation energy for GBS-disl creep
   params.gbs.V_lt1250 = 15e-6  ; % activation volume
   params.gbs.p_lt1250 = 2 ; % grain size exponent
   params.gbs.n_lt1250 = 3.5 ; % stress exponent
   params.gbs.r_lt1250 = 0 ; % water fugacity exponent
   params.gbs.alf_lt1250 = 35 ; % melt factor
   params.gbs.phi_c_lt1250 = phi_c(3);
   params.gbs.x_phi_c_lt1250 = x_phi_c(3);
   
   params.gbs.A_gt1250 = 4.7e10  ; % preexponential for GBS-disl creep
   params.gbs.Q_gt1250 = 600e3  ; % activation energy for GBS-disl creep
   params.gbs.V_gt1250 = 15e-6  ; % activation volume
   params.gbs.p_gt1250 = 2 ; % grain size exponent
   params.gbs.n_gt1250 = 3.5 ; % stress exponent
   params.gbs.r_gt1250 = 0 ; % water fugacity exponent
   params.gbs.alf_gt1250 = 35 ; % melt factor
   params.gbs.phi_c_gt1250=phi_c(3);
   params.gbs.x_phi_c_gt1250=x_phi_c(3);
   
end

function params = load_LH12_flowlaw_constants(phi_c,x_phi_c)
   
% from Hansen et al, 2012 (LH12)=======================
% Coble Diffusion creep (GB)
  params.diff.A = 10^7.6 ; % preexponential for coble diffusion creep
  params.diff.Q = 375e3 ;% activation energy for coble diffusion creep
  params.diff.V = 10e-6 ; % activation volume for coble diff creep
  params.diff.p = 3 ; % grain size exponent
  params.diff.alf = 25 ; % melt factor
  params.diff.r = 0 ; % water fugacity exponent
  params.diff.n = 1 ; % stress exponent
  params.diff.phi_c=phi_c(1);
  params.diff.x_phi_c=x_phi_c(1);
  
% Dislocation creep 
  params.disl.A = 1.1e5 ; % preexponential 
  params.disl.Q = 530e3 ;% activation energy 
  params.disl.V = 15e-6 ; % activation volume 
  params.disl.n = 3.5 ; % stress exponent
  params.disl.p = 0 ; % grain size exponent
  params.disl.alf = 30 ; % melt factor
  params.disl.r = 0 ; % water fugacity exponent
  params.disl.phi_c=phi_c(2);
  params.disl.x_phi_c=x_phi_c(2);
  
% GBS disl accomodated
  params.gbs.A = 10^4.8  ; % preexponential for GBS-disl creep
  params.gbs.Q = 445e3  ; % activation energy for GBS-disl creep
  params.gbs.V = 15e-6  ; % activation volume
  params.gbs.p = 0.73 ; % grain size exponent
  params.gbs.n = 2.9 ; % stress exponent
  params.gbs.alf = 35 ; % melt factor
  params.gbs.r = 0 ; % water fugacity exponent
  params.gbs.phi_c=phi_c(3);
  params.gbs.x_phi_c=x_phi_c(3);
end