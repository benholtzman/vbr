%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates strain rates and viscosities for input state variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VBR= visc_calc_LH2012(VBR)

%% extract state variables and parameters
   T_K = VBR.in.SV.T_K ; % [K]
   P_Pa = 1e9.*(VBR.in.SV.P_GPa) ; % [GPa] to [Pa]
   sig = VBR.in.SV.sig_MPa; % deviatoric stress [MPa]
   d = VBR.in.SV.dg_um ; % [um]
   phi = VBR.in.SV.phi ;
   fH2O=zeros(size(T_K)); % this is a dry flow law  
   params=VBR.in.viscous.LH2012;
   
%% pressur dependent calculation?    
   P_Pa = P_Pa.*(strcmp(params.P_dep_calc,'yes'));    
   
%% calculate water fugacity    
   
   
%% calculate strain rate [1/s]
   sr_tot = 0;
   possible_mechs={'diff';'disl';'gbs'};
      
   for ip = 1:3
       mech=possible_mechs{ip};        
       if isfield(VBR.in.viscous.LH2012,mech)
%         pull out the flow law parameters              
          FLP=params.(mech); 
          
%         calculate strain rate             
          sr = sr_flow_law(T_K,P_Pa,sig,d,phi,fH2O,FLP);              
          sr_tot=sr_tot+sr; 

%         store it          
          VBR.out.viscous.LH2012.(mech).sr=sr; 
          VBR.out.viscous.LH2012.(mech).eta = sig*1e6./sr; % viscosity
       end
       
   end  
  
%% store total composite strain rate and effective viscosity 
   VBR.out.viscous.LH2012.sr_total=sr_tot; % total strain rate 
   VBR.out.viscous.LH2012.eta_total = sig*1e6./sr_tot ; % total viscosity   
  
end
