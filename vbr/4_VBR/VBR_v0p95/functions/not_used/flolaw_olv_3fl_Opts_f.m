function[sr,sr_tot] = flolaw_olv_3fl_Opts_f(T,P,sig_MPa,d,phi,fH2O,params)  
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
%    params        structure with parameters and flags
%
% output: 
%    sr.           strain rate structure. each column corresponds to a different 
%                  deformation mechanism: 
%      sr(1).sr    diffusion (coble) creep
%      sr(2).sr    dislocation creep
%      sr(3).sr    diffusion accomodated gbs
%
%    sr_tot        composite strain rate (i.e. sr(:,1)+sr(:,2)+sr(:,3))
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
  comp = 0; % composition flag or variable -- could be added. 
  n_mech=params.n_mech;    
  sr_tot = zeros(size(T)); 
  P = P.*1e6.*(strcmp(params.P_dep_calc,'yes')); % convert to Pa, zero if we don't want it
  R = 8.314 ; % gas constant 

% get the flowlaw parameter structure, FP  
  FP = build_parameter_struct(fH2O,params); 
  
% loop over mechanisms, calculate strain rates. 
  sr = struct();
  for im=1:n_mech            
%     calculate melt enhancement
      [enhance] = SR_melt_enhancement(phi,FP(im).alf,...
                                 FP(im).speedfac,FP(im).phi_c);
      % 
       meltsr = (1./FP(im).speedfac).* enhance;
       
%     calculate strain rate
      sr(im).sr = (FP(im).A).*(sig_MPa.^FP(im).n)...
                         .*(d.^(-FP(im).p)) .*meltsr ...
                         .*exp(-(FP(im).Q+P.*FP(im).V)./(R.*T))...
                         .*fH2O.^FP(im).r;

%     composite strain rate         
      sr_tot = sr_tot + sr(im).sr;       
  end
  
end

function FP = build_parameter_struct(fH2O,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Builds structure FP for flow law parameters at each thermodynamic condition
% Input: 
%     fH2O   water fugacity [MPa]
% 
% Output:
%  FP.       structure containing flow law parameters. 
%    .A         pre-exponential constant
%    .n         stress exponent
%    .p         grain size exponent
%    .Q         activation energy
%    .V         activation volume
%    .phi_c     critical melt fraction
%    .speedfac  melt enchancement factor
%
% The size of each structure element is (n_T,n_mechanisms) where n_T is the 
% length of the input thermodynamic arrays and n_mechanisms is the number
% of deformation mechanisms. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize constants 
  nmechs = params.n_mech; % number of mechanisms
  FP = struct(); % the structure for flow law constants (FlowParameters=FP)
  
% loop over fields, choosing wet and dry values  
  fields=fieldnames(params.consts); % 
  for i_mech=1:nmechs
      for ifi = 1:numel(fields(:,1))

%       get fieldname    
        name=char(fields(ifi,:)); % fieldname for dry parameters

%       get wet and dry value for current field name      
        dry = params.consts(i_mech).(name)(1);
        wet = params.consts(i_mech).(name)(2);
        
%       set wet or dry depending on water content
        wet_dry=dry .* (fH2O == 0) + wet .* (fH2O > 0);  
              
%       store it in the matrix
        FP(i_mech).(name) =  wet_dry;       
      end      
  end

end