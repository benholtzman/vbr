function [VBR]=checkInput(VBR)
  VBR.status=1;
  VBR.error_message='';
  % make sure the required methods are set:
  if isfield(VBR.in,'anelastic')
    if sum(strncmp('YT_maxwell',VBR.in.anelastic.methods_list,10)) > 0
      if isfield(VBR.in,'viscous') == 0
        VBR.status=0;
        msg='VBR failure: YT_maxwell requires a viscosity method e.g., ';
        msg=strcat(msg,'VBR.in.viscous.methods_list={''LH2012''};');
        msg=strcat(msg,'YT_maxwell will only use the first specified viscosity method');
        VBR.error_message=msg;
      end
    end
  end

end
