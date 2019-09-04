function [VBR] = VBR_spine(VBR)

%% =====================================================================
%% Check VBR Input
%% =====================================================================
  VBR = checkInput(VBR);
  telapsed=struct();
  if VBR.status==0
     disp(VBR.error_message)
     return
  end

%% =====================================================================
%% ELASTIC properties ==================================================
%% =====================================================================
if isfield(VBR.in,'elastic')
  [VBR,telapsed.elastic]=spineGeneralized(VBR,'elastic');
end

%% =====================================================================
%% VISCOUS properties ==================================================
%% =====================================================================

if isfield(VBR.in,'viscous')
  if isfield(VBR.in.viscous,'methods_list')
   [VBR,telapsed.visc]=spineGeneralized(VBR,'viscous');
  end
end

%% =====================================================================
%% ATTENUATION ! =======================================================
%% =====================================================================

if isfield(VBR.in,'anelastic')
   [VBR,telapsed.anelastic]=spineGeneralized(VBR,'anelastic');
end

%% ========================================================================
   VBR.out.computation_time=telapsed; % store elapsed time for each

end
