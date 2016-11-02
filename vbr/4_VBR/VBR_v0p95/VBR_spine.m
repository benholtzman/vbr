function[VBR] = VBR_spine(VBR) 
%% =====================================================================
%% ELASTIC properties ==================================================
%% =====================================================================

if isfield(VBR.in,'elastic') 
   telapsed.elastic=tic;
   
   methods_list=VBR.in.elastic.methods_list; % list of methods to use
   
%% anharmonic, unrelaxed calculation   
   if sum(strncmp('anharmonic',methods_list,10)) > 0 || ...
      sum(strncmp('poro_Takei',methods_list,10)) > 0  
  
      if isfield(VBR.in.elastic,'anharmonic')==0   
           VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); 
      end
      
     [VBR] = el_calc_Gu_0(VBR); % get reference shear modulus
     [VBR] = el_ModUnrlx_dTdP_f(VBR) ; % project to T, P of interest      
   end
   
%% poroelastic effect of melt   
   if sum(strncmp('poro_Takei',methods_list,10)) > 0
       if isfield(VBR.in.elastic,'poro_Takei')==0
           VBR.in.elastic.poro_Takei=Params_Elastic('poro_Takei');
       end
      [VBR] = el_ModUnrlx_MELT_f(VBR) ; 
   end
      
%% stixrude and lithgow-bertelloni scaling   
   if sum(strncmp('SLB2005',methods_list,16)) > 0    
       [VBR] = el_Vs_SnLG_f(VBR); % Stixrude and Lithgow-Bertelloni
   end
   
   telapsed.elastic=toc(telapsed.elastic);
end

%% =====================================================================
%% VISCOUS properties ==================================================
%% =====================================================================

if isfield(VBR.in,'viscous')
   telapsed.visc=tic;
   
   methods_list=VBR.in.viscous.methods_list; % list of methods to use
   
   if sum(strncmp('HK2003',methods_list,6)) > 0 
      if isfield(VBR.in.viscous,'HK2003')==0
          VBR.in.viscous.HK2003=Params_Viscous('HK2003');
      end 
      VBR = visc_calc_HK2003(VBR);       
   end
   if sum(strncmp('LH2012',methods_list,6)) > 0       
      if isfield(VBR.in.viscous,'LH2012')==0
          VBR.in.viscous.LH2012=Params_Viscous('LH2012');
      end 
      VBR = visc_calc_LH2012(VBR);       
   end
   
   
   telapsed.visc=toc(telapsed.visc);
end

%% =====================================================================
%% ATTENUATION ! =======================================================
%% =====================================================================

if isfield(VBR.in,'anelastic')

   methods_list=VBR.in.anelastic.methods_list; % list of methods to use

%  Extended Burgers method 
   if sum(strncmp('eBurgers',methods_list,8)) > 0
      if isfield(VBR.in.anelastic,'eBurgers')==0
          VBR.in.anelastic.eBurgers=Params_Anelastic('eBurgers');
      end
      if strcmp(VBR.in.anelastic.eBurgers.method,'PointWise') 
          telapsed.eBurgers=tic;
          VBR=Q_eBurgers_f(VBR);
          telapsed.eBurgers=toc(telapsed.eBurgers);
      elseif strcmp(VBR.in.anelastic.eBurgers.method,'FastBurger')
          telapsed.eBurgers=tic;
          [VBR]=Q_eFastBurgers(VBR) ; 
          telapsed.eBurgers=toc(telapsed.eBurgers);   
      end
   end

%  Andrade Pesudo-Period method
   if sum(strncmp('AndradePsP',methods_list,10)) > 0
      telapsed.AndradePsP=tic;
      if isfield(VBR.in.anelastic,'AndradePsP')==0
         VBR.in.anelastic.AndradePsP=Params_Anelastic('AndradePsP');
      end
      [VBR]=Q_Andrade_PseudoP_f(VBR) ; 
      telapsed.AndradePsP=toc(telapsed.AndradePsP);
   end
   
end

%% ========================================================================
   VBR.out.computation_time=telapsed; % store elapsed time for each

end
    
