function [VBR,telapsed]=spineElastic(VBR)
  % calculates user-selected elastic methods

  telapsed.elastic=tic;

  methods_list=VBR.in.elastic.methods_list; % list of methods to use

  %% anharmonic, unrelaxed calculation
  if sum(strncmp('anharmonic',methods_list,10)) > 0 || ...
     sum(strncmp('poro_Takei',methods_list,10)) > 0
     % if you are not calculating the anharmonic properties,
     % still use the anharmonic reference moduli:
     if isfield(VBR.in.elastic,'anharmonic')==0
          VBR.in.elastic.anharmonic=Params_Elastic('anharmonic');
     end

    [VBR] = el_calc_Gu_0(VBR); % get reference shear modulus from the anharmonic method
    [VBR] = el_ModUnrlx_dTdP_f(VBR) ; % calculate anharmonic scaling with T, P of interest
  end

  %% poroelastic effect of melt
  if sum(strncmp('poro_Takei',methods_list,10)) > 0
      if isfield(VBR.in.elastic,'poro_Takei')==0
          VBR.in.elastic.poro_Takei=Params_Elastic('poro_Takei');
      end
     [VBR] = el_ModUnrlx_MELT_f(VBR) ;
  end

  %% stixrude and lithgow-bertelloni scaling
  if sum(strncmp('SLB2005',methods_list,7)) > 0
      [VBR] = el_Vs_SnLG_f(VBR); % Stixrude and Lithgow-Bertelloni
  end

  telapsed.elastic=toc(telapsed.elastic);
end
