function [VBR,telapsed]=spineAnelastic_proposed(VBR)
  % calculates user-selected anelastic methods

  methods_list=VBR.in.anelastic.methods_list; % list of methods to use
  empty_params=Params_Anelastic('');
  possible_methods=empty_params.possible_methods; % list of allowed methods
  telapsed=struct();

  for i_method = 1:numel(methods_list)
    meth=methods_list{i_method};
    if any(strcmp(possible_methods,meth))      
      % load parameters for this method if not set:
      if isfield(VBR.in.anelastic,meth)==0
          VBR.in.anelastic.(meth)=Params_Anelastic(meth);
      end

      % call the function specificed by this method
      func_name=VBR.in.anelastic.(meth).func_name;
      telapsed.(meth)=tic;
      [VBR] = feval(func_name,VBR);
      telapsed.(meth)=toc(telapsed.(meth));
    end
  end

end
