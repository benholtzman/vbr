function [VBR] = Q_eBurgers_decider(VBR)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % [VBR] = Q_eBurgers(VBR)
  %
  % wrapper for eBurger methods. decides which eBurgers method to call: either
  % PointWise or FastBurger method. FastBurger only works with the high temp
  % background, use PointWise if using a peak.
  %
  % Parameters:
  % ----------
  % VBR   the VBR structure
  %
  % Output:
  % -----
  % VBR   the VBR structure, with VBR.out.anelastic.eBurgers structure
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if strcmp(VBR.in.anelastic.eBurgers.method,'PointWise')
    VBR=Q_eBurgers_f(VBR);
  elseif strcmp(VBR.in.anelastic.eBurgers.method,'FastBurger')
    [VBR]=Q_eFastBurgers(VBR) ;
  end
end
