%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates strain rates and viscosities for input state variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VBR= visc_calc_HK2003(VBR)

%% extract state variables and parameters
   T_K = VBR.in.SV.T_K ; % [K]
   P_Pa = 1e9.*(VBR.in.SV.P_GPa) ; % [GPa] to [Pa]
   sig = VBR.in.SV.sig_MPa; % deviatoric stress [MPa]
   d = VBR.in.SV.dg_um ; % [um]
   phi = VBR.in.SV.phi ;
   Ch2o=VBR.in.SV.Ch2o; % in [PPM]!
   ch2o=VBR.in.viscous.HK2003.ch2o_o; % effectively dry below this [PPM]
   params=VBR.in.viscous.HK2003;
   
%% pressur dependent calculation?    
   P_Pa = P_Pa.*(strcmp(params.P_dep_calc,'yes'));    
   
%% calculate water fugacity    
   fH2O=water_fugacity(Ch2o,ch2o,P_Pa,T_K); % [MPa]
   VBR.in.SV.Fh2o=fH2O;
  
%% calculate strain rate [1/s]
   sr_tot = 0;
   possible_mechs={'diff';'disl';'gbs'};
      
   for ip = 1:3
       mech=possible_mechs{ip};        
       if isfield(VBR.in.viscous.HK2003,mech)
%         prep the flow law parameters              
          FLP=prep_constants(fH2O,T_K,params.(mech),mech); 
          
%         calculate strain rate   
          sr = sr_flow_law(T_K,P_Pa,sig,d,phi,fH2O,FLP);    
          sr_tot=sr_tot+sr; 
          
          VBR.out.viscous.HK2003.(mech).sr=sr; 
          VBR.out.viscous.HK2003.(mech).eta = sig*1e6./sr; % viscosity
       end
       
   end  
  
%% store total composite strain rate and effective viscosity 
   VBR.out.viscous.HK2003.sr_total=sr_tot; % total strain rate 
   VBR.out.viscous.HK2003.eta_total = sig*1e6./sr_tot ; % total viscosity   
  
end

function FLP=prep_constants(fH2O,T_K,params,mech)
   FLP = struct(); 
   fields={'A';'n';'p';'r';'Q';'V';'phi_c';'alf';'x_phi_c'}; 
   for ifi = 1:numel(fields(:,1))
%      get fieldname
       name=char(fields(ifi,:)); 

       if strcmp(mech,'diff') || strcmp(mech,'disl')          
          name_wet = [name '_wet'];          
          dry = params.(name);
          wet = params.(name_wet);
%         set wet or dry depending on water content
          wet_dry=dry .* (fH2O == 0) + wet .* (fH2O > 0);   
%         store it in the matrix
          FLP.(name) =  wet_dry;
          
       elseif strcmp(mech,'gbs')
           
          name_gt = [name '_gt1250']; 
          name_lt = [name '_lt1250'];           
          gt1250 = params.(name_gt);          
          lt1250 = params.(name_lt);          
%         set wet or dry depending on water content
          val = gt1250 .* (T_K-273 >= 1250) + lt1250 .* (T_K-273<1250);   
%         store it in the matrix
          FLP.(name) =  val;           
       end
                     
   end
end
          
