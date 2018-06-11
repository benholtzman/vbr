%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [VBR]=Q_eFastBurgers(VBR) 
% c. havlin, 04-2016
%
% calculates material properties for extended burgers model using the
% FastBurger integration algorithm. Rather than integrate over relaxation 
% period for every thermodynamic state (as in Q_eBurgers_f.m), the 
% FastBurger algorithm first integrates over the entire range of relaxation 
% periods of all input thermodynamic states, then extracts the relevant 
% range for each thermodynamic state. See Notes for detailed description.
%  
% input:
%      VBR.   VBR structure with state variables and eBurger settings
%
% output:
%      VBR.eBurgers.    
%                  .J1     J1 comliance [1/Pa]
%                  .J2     J2 compliance [1/Pa]
%                  .Q      (J1/J2)
%                  .Qinv   attenuation (Q^{-1})
%                  .M      modulus [Pa]
%                  .V      relaxed soundspeed [m/s]
%                  .Vave   averaged V over all frequencies [m/s]
%
% A note on dimensions: arrays for frequency dependent variables will have 
% an extra dimension for each frequency supplied. i.e., V(4,5,2) will be
% the velocity at thermodynamic point (4,5) and frequency freq_vec(2). Vave
% is the only output variable that is not frequency dependent, thus has the
% same size as incoming thermodynamic state variables. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[VBR]=Q_eFastBurgers(VBR) 
%% ===========================
%% read in thermodynamic state
%% ===========================
   f_vec = VBR.in.SV.f ;   
   phi =  VBR.in.SV.phi ; 
   if isfield(VBR.in.elastic,'poro_Takei')
     Mu = VBR.out.elastic.poro_Takei.Gu ;   
   else if isfield(VBR.in.elastic,'anharmonic')
     Mu = VBR.out.elastic.anharmonic.Gu ;    
   end    
   T_K_mat = VBR.in.SV.T_K ;
   P_Pa_mat = VBR.in.SV.P_GPa.*1e9 ; % convert pressure GPa to Pa = GPa*1e9
   rho_mat = VBR.in.SV.rho ;
   d_mat = VBR.in.SV.dg_um ; % microns grain size
   
   w_vec = 2*pi.*f_vec ; % period
   Ju_mat = 1./Mu ; % unrelaxed compliance
   
%% =================
%% allocate matrices  
%% =================
%  frequency is added as a new dimension at end of array.
%  e.g., if size(T_K_mat)=n1,n2 then the frequency dependent variables will
%        have size(Q)=n1,n2,nfreq. 
   nfreq = numel(f_vec); 
   Jz=zeros(size(Ju_mat)); 
   J1 = proc_add_freq_indeces(Jz,nfreq); 
   J2 = J1; Q = J1; Qinv = J1; M = J1; V = J1;             
   Vave = Jz; 

%% read in reference values 
   Burger_params=VBR.in.anelastic.eBurgers;
   TR = Burger_params.TR ;% Kelvins
   PR = Burger_params.PR *1e9; % GPa to Pa
   dR = Burger_params.dR ; % microns grain size  
   E = Burger_params.E ; % J/mol
   R = Burger_params.R ;
   Vstar = Burger_params.Vstar ; % m^3/mol (Activation Volume? or molar volume?)
   m = Burger_params.m ;
%  Q settings (Jackson n Faul 2010, table 1) : 
   alf = Burger_params.alf ; % is this the same as n in Andrade ? a constant 
   Delta = Burger_params.Delta ;%1.4 ;% ; % relaxation strength..   
   Tau_LR = Burger_params.Tau_LR ;
   Tau_HR = Burger_params.Tau_HR ;
   
%  melt effects
   alpha = Burger_params.melt_alpha ; 
   phi_c = Burger_params.phi_c ;
   x_phi_c = Burger_params.x_phi_c ; 

%% scaling 
   scale_mat = ((d_mat./dR).^m).*exp((E/R).*(1./T_K_mat-1/TR)) ...
                            .*exp((Vstar/R).*(P_Pa_mat./T_K_mat-PR/TR)) ; 
    % more like viscosity, so divide by melt factor 
    % (Xtilde is like strain rate!, where you multiply by rate factor) 

% EFFECT OF MELT, hypothetical: 
% sharper response than creep at phi_c, but less sensitive at higher phi (for HTB only) 
scale_mat = scale_mat.*x_phi_c ; % to account for lack of truly melt free samples and the drop at the onset of melting.  
[scale_mat_prime] = sr_melt_enhancement(phi,alpha,x_phi_c,phi_c) ;
scale_mat = scale_mat./scale_mat_prime ; 


% define global Tau range to integrate over
  Tau_L = Tau_LR.*scale_mat ; 
  Tau_H = Tau_HR.*scale_mat ; 
  
  nglobvec0=Burger_params.nTauGlob;
  minTau_L=min(Tau_L(1:numel(Tau_L)));
  maxTau_H=max(Tau_H(1:numel(Tau_H)));
  nTau = numel(Tau_L);   
  Tau_glob_vec0(1:nglobvec0) = logspace(log10(minTau_L),log10(maxTau_H),nglobvec0);

% make sure each Tau_L, Tau_H appears in vector exactly  
  for i_th = 1:nTau; 
       Tau_L0 = Tau_L(i_th);
       Tau_H0 = Tau_H(i_th);
       
       if isempty(find(Tau_glob_vec0==Tau_L0));
           Tau_glob_vec0 = [Tau_glob_vec0 Tau_L0];
       end
       if isempty(find(Tau_glob_vec0==Tau_H0));
           Tau_glob_vec0 = [Tau_glob_vec0 Tau_H0];
       end
  end

% sort the new vector, track indeces  
  [Tau_glob_vec,TauIndx]=sort(Tau_glob_vec0);  
  nglobvec=numel(Tau_glob_vec); 
  
% integrate the portions that depend only on frequency, cumulatively
  ints=repmat(struct('J1', nglobvec,'J2', nglobvec), nfreq, 1 );      
  for iw = 1:nfreq
   w = w_vec(iw); % for now   
   
   ints(iw).J1= cumtrapz(Tau_glob_vec,Tau_glob_vec.^(alf-1)./(1+w^2.*Tau_glob_vec.^2)) ;   
   ints(iw).J2 = cumtrapz(Tau_glob_vec,Tau_glob_vec.^alf./(1+w^2.*Tau_glob_vec.^2)) ;
  end
  
% now go through each i_th, find the interval of interest for each and 
% calculate the modulus 
for i_th = 1:nTau
%   get current i_th values    
    scale = scale_mat(i_th);
    Ju = Ju_mat(i_th) ;
    rho = rho_mat(i_th) ;
    
%   the bounds to find            
    Tau_L0 = Tau_L(i_th);
    Tau_H0 = Tau_H(i_th);
          
%   find the bounds
    iLow = find(Tau_glob_vec==Tau_L0); 
    iHigh = find(Tau_glob_vec==Tau_H0);
       
%   loop over frequency
    for iw=1:nfreq        
        i_glob = i_th + (iw - 1) * nTau; % the linear index of the arrays with
                                         % a frequency index
                
        w = w_vec(iw) ;
        Tau_MR = 10^5.2 ;
        Tau_M = Tau_MR.*scale;
        
        int_J1 = ints(iw).J1;
        d_integral=alf*(int_J1(iHigh)-int_J1(iLow))./(Tau_H0^alf - Tau_L0^alf);
        J_int_1_0 = (1+Delta*d_integral) ;
        
        int_J2 = ints(iw).J2;
        d_integral=alf*(int_J2(iHigh)-int_J2(iLow))./(Tau_H0^alf - Tau_L0^alf);
        J_int_2_0 = (w*Delta*d_integral + 1/(w*Tau_M)) ;
        
        J1(i_glob)=Ju * J_int_1_0;
        J2(i_glob)=Ju * J_int_2_0;
        
        Q(i_glob) = J1(i_glob)./J2(i_glob) ;
        Qinv(i_glob) = 1./Q(i_glob) ;
        M(i_glob) = (J1(i_glob).^2 + J2(i_glob).^2).^(-0.5) ;
        V(i_glob) = sqrt(M(i_glob)./rho) ;
        
        Vave(i_th) = Vave(i_th) + V(i_glob);
    end
       
end

%% WRITE VBR
   VBR.out.anelastic.eBurgers.J1 = J1; 
   VBR.out.anelastic.eBurgers.J2 = J2; 
   VBR.out.anelastic.eBurgers.Q = Q; 
   VBR.out.anelastic.eBurgers.Qinv = Qinv; 
   VBR.out.anelastic.eBurgers.M=M; 
   VBR.out.anelastic.eBurgers.V=V; 
   VBR.out.anelastic.eBurgers.Vave = Vave./nfreq;

end
    
