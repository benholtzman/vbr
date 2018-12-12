function [VBR]=checkInput(VBR)
  % checkInput() looks for some dependency errors. e.g., the YT_maxwell method
  % requires a viscous method to run, all the anelastic calculation require a
  % pure elastic calculation.

  % initialize error storage
  VBR.status=1;
  VBR.error_message='';

  % build requirements lists
  Reqs={'anelastic','YT_maxwell','viscous';
        'anelastic','YT_maxwell','elastic';
        'anelastic','eBurgers','elastic';
        'anelastic','AndradePsP','elastic';
        'anelastic','YT2016_solidus','elastic';
        'anelastic','YT2016_solidus','SV.Tsolidus_K'};
        
  % list of messages to display for each of those requirements
  Msgs={'YT_maxwell requires a viscosity method. e.g., VBR.in.viscous.methods_list={''LH2012''};';
        'YT_maxwel requires an elasticity method. e.g., VBR.in.elastic.methods_list={''anharmonic''}';
        'eBurgers requires an elasticity method. e.g., VBR.in.elastic.methods_list={''anharmonic''}';
        'AndradePsP requires an elasticity method. e.g., VBR.in.elastic.methods_list={''anharmonic''}';
        'YT2016_solidus requires an elasticity method. e.g., VBR.in.elastic.methods_list={''anharmonic''}';
        'YT2016_solidus requires a solidus state variable e.g., VBR.in.SV.Tsolidus_K=1200*ones(size(T))'};

  % loop over requirements, check them.
  for ri = 1:size(Reqs)(1)
      typ=Reqs{ri,1}; % general method field
      if isfield(VBR.in,typ)
        % check if the method is invoked
        method=Reqs{ri,2};
        methods=VBR.in.(typ).methods_list;
        if sum(strncmp(method,methods,8)) > 0
          % check the required fields/methods for this method
          fieldvars=strsplit(Reqs{ri,3},'.');
          if isfield(VBR.in,fieldvars{1})==0
            VBR.status=0;
          else
            if numel(fieldvars)>1
              if isfield(VBR.in.(fieldvars{1}),fieldvars{2})==0
                VBR.status=0;
              end
            else
          end

          % save the message & exit if error caught
          if VBR.status==0
            msg=Msgs{ri};
            msg=strcat("\n!!! VBR FAILURE !!!\n",msg);
            msg=strcat(msg,"\n!!! VBR FAILURE !!!\n");
            VBR.error_message=msg;
            return;
          end
        end
      end

  end

end
