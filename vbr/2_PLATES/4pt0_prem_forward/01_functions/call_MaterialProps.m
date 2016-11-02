%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Vark = call_MaterialProps(Vark,Rho_0,Kc_0,Cp_0,settings,z,dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                call_MaterialProps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to do some trickery... better to calculate material properties only
% on interior, nonghosts. the ghosts mess up the pressure integration. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zero the things
  nz = numel(Vark.T);
  Vark.rho = zeros(nz,1);
  Vark.cp = zeros(nz,1);
  Vark.kc = zeros(nz,1);
  Vark.P = zeros(nz,1);

% want to use the staggered values -- the integration for P is easier
% using cell edges.  
  Ts = stag(Vark.T)+273; rho = stag(Rho_0);
  kc = stag(Kc_0); cp = stag(Cp_0); 
  zs = stag(z);   
% calculate material props     
  [rho,cp,kc,P_hyd]=MaterialProperties(rho,kc,cp,Ts,zs,settings.P0,...
                     settings.dTdz_ad,settings.Flags.PropType); 
                 
% fill the ghosts   
  Vark.cp = addghosts(stag(cp));
  Vark.Kc = addghosts(stag(kc));
  Vark.rho = addghosts(stag(rho));
  Vark.P = addghosts(stag(P_hyd)); 
  
  if isreal(Vark.rho)==0 || isreal(Vark.cp)==0 || isreal(Vark.Kc)==0
      disp('   wtf')
  end
   
end