%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [VBR] = Q_eBurgers(VBR)
% wrapper for eBurger methods: decides which eBurgers method to call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [VBR] = Q_eBurgers_decider(VBR)
  if strcmp(VBR.in.anelastic.eBurgers.method,'PointWise')
    VBR=Q_eBurgers_f(VBR);
  elseif strcmp(VBR.in.anelastic.eBurgers.method,'FastBurger')
    [VBR]=Q_eFastBurgers(VBR) ;
  end
end
