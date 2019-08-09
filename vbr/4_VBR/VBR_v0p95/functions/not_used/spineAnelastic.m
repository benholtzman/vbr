function [VBR,telapsed]=spineAnelastic(VBR)
  % calculates user-selected anelastic methods

  methods_list=VBR.in.anelastic.methods_list; % list of methods to use
  empty_params=Params_Anelastic('');
  possible_methods=empty_params.possible_methods; % list of allowed methods
  telapsed=struct();

  for i_method = 1:numel(methods_list)
    meth=methods_list{i_method};
    if any(strcmp(possible_methods,meth))
      % load parameters for this method, keep parameters set by user
      meth_params=Params_Anelastic(meth);

      if isfield(VBR.in.anelastic,meth)
        % user set parameters, check which they set
        fldz=fieldnames(meth_params);
        for ifield=1:numel(fldz)
          if ~isfield(VBR.in.anelastic.(meth),fldz{ifield})
            % user didn't set this parameter, take the default
            VBR.in.anelastic.(meth).(fldz{ifield})=meth_params.(fldz{ifield});
          end
        end
      else
        % user did not set any parameters, use all defaults
        VBR.in.anelastic.(meth)=meth_params;
      end

      % call the function specificed by this method
      func_name=VBR.in.anelastic.(meth).func_name;
      telapsed.(meth)=tic;
      [VBR] = feval(func_name,VBR);
      telapsed.(meth)=toc(telapsed.(meth));
    end
  end

end
