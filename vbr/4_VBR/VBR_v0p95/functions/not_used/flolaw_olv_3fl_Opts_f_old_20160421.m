function[sr,sr_tot] = flolaw_olv_3fl_Opts_f(T,P,sig_MPa,d,phi,fH2O,FlowLawOpts)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates strain rates at the given thermodynamic conditions
%
% input:
%    T             temperature in Kelvins
%    P             pressure in MPa (converted to Pa below)
%    sig           stress in MPa
%    d             grain size in microns
%    phi           melt fraction (i.e. 0<phi<1)
%    fH2O          water fugacity in MPa
%    FlowLawOpts   structure with parameters and flags
%
% output: 
%    sr            strain rate matrix. each column corresponds to a different 
%                  deformation mechanism: 
%      sr(:,1)     diffusion (coble) creep
%      sr(:,2)     dislocation creep
%      sr(:,3)     diffusion accomodated gbs
%
%    sr_tot        composite strain rate (i.e. sr(:,1)+sr(:,2)+sr(:,3))
% 
% C. Havlin 10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
  comp = 0; % composition flag or variable -- could be added. 
  n_mech=FlowLawOpts.n_mech; 
  [nTS1,nTS2,nTS3] = size(T); % number of thermodynamic states
  sr = zeros(nTS1,nTS2,nTS3,n_mech); % strain rate matrix (column for each mechanism) 
  sr_tot = zeros(nTS1,nTS2,nTS3); 
  P = P.*1e6.*(strcmp(FlowLawOpts.P_dep_calc,'yes')); % convert to Pa, zero if we don't want it
  R = 8.314 ; % gas constant 

% get the flowlaw parameter structure, FP  
  FP = initialize_flow_laws(T,fH2O,nTS1,nTS2,nTS3,comp,FlowLawOpts); 
  
% loop over mechanisms, calculate strain rates. 
  sr = struct();
  for im=1:n_mech            
%     calculate melt enhancement
      [enhance] = SR_melt_enhancement(phi,FP.alf(:,:,:,im),...
                                 FP.speedfac(im),FP.phi_c(im));
      % 
       meltsr = (1./FP.speedfac(im)).* enhance;
       
%     calculate strain rate
      sr(im).sr(:,:,:) = (FP.A(:,:,:,im)).*(sig_MPa.^FP.n(:,:,:,im))...
                         .*(d.^(-FP.p(:,:,:,im))) .*meltsr ...
                         .*exp(-(FP.Q(:,:,:,im)+P.*FP.V(:,:,:,im))./(R.*T))...
                         .*fH2O.^FP.r(:,:,:,im);

%     composite strain rate         
      sr_tot = sr_tot + sr(im).sr(:,:,:); 
      sr(im).sr= squeeze(sr(im).sr);
  end
  
  
  sr_tot = squeeze(sr_tot);

end


function FP = initialize_flow_laws(T,fH2O,nTS1,nTS2,nTS3,comp,Params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Builds structure FP for flow law parameters at each thermodynamic condition
% Input: 
%     T      temperature [K]
%     fH2O   water fugacity [MPa]
%     phi    melt fraction
%     comp   composition (not used right now...)
% 
% Output:
%  FP.       structure containing flow law parameters
%    .A      pre-exponential constant
%    .n      stress exponent
%    .p      grain size exponent
%    .Q      activation energy
%    .V      activation volume
%    .melt   melt enhancement
%
% The size of each structure element is (n_T,n_mechanisms) where n_T is the 
% length of the input thermodynamic arrays and n_mechanisms is the number
% of deformation mechanisms. 
% 
% as of now, the columns of each element correspond ot the following
% deformation mechanisms:
%    1 = coble (diffusion), 2 = gbs, 3 = dislocation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parameters
  StudyNameDry=Params.StudyNameDry;
  StudyNameWet=Params.StudyNameWet;  
  phi_c = Params.phi_c ;  
  speedfacdiff = Params.speedfac;
  
% load all flow law constants
  Rhe = flow_law_consts; 

% initialize constants 
  nmechs = Params.n_mech; % number of mechanisms
  FP.A = (zeros(nTS1,nTS2,nTS3,nmechs)); % pre-exponential
  FP.n = FP.A ; % stress exponent
  FP.p = FP.A ; % grain size exponent
  FP.r = FP.A; % water fugacity exponent
  FP.Q = FP.A; % activation energy
  FP.V = FP.A; % activation volume
  FP.alf = FP.A; % standard melt fraction exponential factor
    
  fields = fieldnames(FP); % the above field names
 
  
% select study for GBS  
  GBSStudyNameDry='LH12'; % only use LH12 for GBS -- HK03 has an annoying step function
   % also, gbs under wet conditions is unknown. currently, we only use dry 
   % gbs, even when using wet conditions for the other mechanisms. 
   
   
  for ifi=1:numel(fields(:,1)) % loop over field names
      for i_mech = 1:nmechs % loop over deformation mechanisms

%       get fieldnames        
        name=char(fields(ifi,:)); % fieldname for dry parameters
        wet_name=[name '_wet']; % fieldname for wet parameters
        
%       extract the constants        
        dry = [Rhe.(['olv_' StudyNameDry '_diff']).(name) ...
              Rhe.(['olv_' StudyNameDry '_disl']).(name) ...
              Rhe.(['olv_' GBSStudyNameDry '_gbs']).(name)];              
        wet = [Rhe.(['olv_' StudyNameWet '_diff']).(wet_name) ...
               Rhe.(['olv_' StudyNameWet '_disl']).(wet_name) ...
               Rhe.(['olv_' GBSStudyNameDry '_gbs']).(name)]; 
              
        
%      set wet or dry dependeing on wetness
       wet_dry=dry(i_mech) .* (fH2O == 0) + wet(i_mech) .* (fH2O > 0);  
              
       FP.(name)(:,:,:,i_mech)=  wet_dry;       
      end      
  end
  
%%%%%%%%%%%%%%%%%%%
% EFFECTS of MELT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assumes: GBS is less sensitive to melt than diff., but still a bit. 
%          no sensitivity of disl. to phic.

  FP.speedfac = [speedfacdiff; 1.0; speedfacdiff/2];
  FP.phi_c = [phi_c; phi_c; phi_c]; % look!! a row of phi_cs !     

end


function Rhe = flow_law_consts
% =========================================================================
% loads constants for all deformation mechanisms and conditions into the 
% Rhe structure. 
%
% Currently uses Olivine Fo90 flow law parameters from:
%           Hirth and Kohlstedt, 2003 (HK03)
%           Hansen et al, 2012 (LH12) 
% ========================================================================= 

% Coble Diffusion creep (GB)
%  dry
   Rhe.olv_HK03_diff.A = 1.5e9 ; % preexponential for coble diffusion creep
   Rhe.olv_HK03_diff.Q = 375e3 ;% activation energy for coble diffusion creep
   Rhe.olv_HK03_diff.V = 10e-6 ; % activation volume for coble diff creep 
   Rhe.olv_HK03_diff.p = 3 ; % grain size exponent
   Rhe.olv_HK03_diff.n = 1 ; % stress exponent
   Rhe.olv_HK03_diff.alf = 25 ; % melt factor
   Rhe.olv_HK03_diff.r = 0 ; % water fugacity exponent
%  wet   
   Rhe.olv_HK03_diff.A_wet = 2.5e7 ; % preexponential for coble diffusion creep
   Rhe.olv_HK03_diff.Q_wet = 375e3 ;% activation energy for coble diffusion creep
   Rhe.olv_HK03_diff.V_wet = 10e-6 ; % activation volume for coble diff creep
   Rhe.olv_HK03_diff.p_wet = 3 ; % grain size exponent
   Rhe.olv_HK03_diff.n_wet = 1 ; % stress exponent
   Rhe.olv_HK03_diff.alf_wet = 25 ; % melt factor
   Rhe.olv_HK03_diff.r_wet = 0.7 ; % water fugacity exponent
  
% Dislocation creep 
  Rhe.olv_HK03_disl.A = 1.1e5 ; % preexponential  
  Rhe.olv_HK03_disl.Q = 530e3 ;% activation energy 
  Rhe.olv_HK03_disl.V = 15e-6 ; % activation volume (HK03 doesn't report, using value from LH12 here)
  Rhe.olv_HK03_disl.n = 3.5 ; % stress exponent
  Rhe.olv_HK03_disl.p = 0; % grain size exponent
  Rhe.olv_HK03_disl.alf = 30 ; % melt factor
  Rhe.olv_HK03_disl.r = 0 ; % water fugacity exponent
  Rhe.olv_HK03_disl.A_wet = 1600.0 ; % preexponential  
  Rhe.olv_HK03_disl.Q_wet = 520e3 ;% activation energy 
  Rhe.olv_HK03_disl.V_wet = 22e-6 ; % activation volume 
  Rhe.olv_HK03_disl.n_wet = 3.5 ; % stress exponent
  Rhe.olv_HK03_disl.p_wet = 0; % grain size exponent
  Rhe.olv_HK03_disl.alf_wet = 30 ; % melt factor
  Rhe.olv_HK03_disl.r_wet = 1.2 ; % water fugacity exponent
  
% GBS disl accomodated (THIS ONE WILL BREAK RIGHT NOW!)
  Rhe.olv_HK03_gbs.A_lt1250 = 6500  ; % preexponential for GBS-disl creep
  Rhe.olv_HK03_gbs.Q_lt1250 = 400e3  ; % activation energy for GBS-disl creep
  Rhe.olv_HK03_gbs.V_lt1250 = 15e-6  ; % activation volume
  Rhe.olv_HK03_gbs.A_gt1250 = 4.7e10  ; % preexponential for GBS-disl creep
  Rhe.olv_HK03_gbs.Q_gt1250 = 600e3  ; % activation energy for GBS-disl creep
  Rhe.olv_HK03_gbs.V_gt1250 = 15e-6  ; % activation volume
  Rhe.olv_HK03_gbs.p = 2 ; % grain size exponent
  Rhe.olv_HK03_gbs.n = 2 ; % stress exponent
  Rhe.olv_HK03_gbs.alf = 35 ; % melt factor
  Rhe.olv_HK03_gbs.r = 0 ; % water fugacity exponent
  
% from Hansen et al, 2012 (LH12)=======================
% Coble Diffusion creep (GB)
  Rhe.olv_LH12_diff.A = 10^7.6 ; % preexponential for coble diffusion creep
  Rhe.olv_LH12_diff.Q = 375e3 ;% activation energy for coble diffusion creep
  Rhe.olv_LH12_diff.V = 10e-6 ; % activation volume for coble diff creep
  Rhe.olv_LH12_diff.p = 3 ; % grain size exponent
  Rhe.olv_LH12_diff.alf = 25 ; % melt factor
  Rhe.olv_LH12_diff.r = 0 ; % water fugacity exponent
  Rhe.olv_LH12_diff.n = 1 ; % stress exponent
  
% Dislocation creep 
  Rhe.olv_LH12_disl.A = 1.1e5 ; % preexponential 
  Rhe.olv_LH12_disl.Q = 530e3 ;% activation energy 
  Rhe.olv_LH12_disl.V = 15e-6 ; % activation volume 
  Rhe.olv_LH12_disl.n = 3.5 ; % stress exponent
  Rhe.olv_LH12_disl.p = 0 ; % grain size exponent
  Rhe.olv_LH12_disl.alf = 30 ; % melt factor
  Rhe.olv_LH12_disl.r = 0 ; % water fugacity exponent
  
% GBS disl accomodated
  Rhe.olv_LH12_gbs.A = 10^4.8  ; % preexponential for GBS-disl creep
  Rhe.olv_LH12_gbs.Q = 445e3  ; % activation energy for GBS-disl creep
  Rhe.olv_LH12_gbs.V = 15e-6  ; % activation volume
  Rhe.olv_LH12_gbs.p = 0.73 ; % grain size exponent
  Rhe.olv_LH12_gbs.n = 2.9 ; % stress exponent
  Rhe.olv_LH12_gbs.alf = 35 ; % melt factor
  Rhe.olv_LH12_gbs.r = 0 ; % water fugacity exponent
end
