function plot_Bayes_surf(P_mod)
  % plots 3d slices of Bayesian Posterior Distribution
  figure()
  Y=P_mod.States.Tpot;
  X=P_mod.States.phi;
  Z=P_mod.States.gs/1e6;
  V=P_mod.Posterior_scaled;

  i_best=find(P_mod.Posterior_scaled(:)==max(P_mod.Posterior_scaled(:)));
  xslice = [X(i_best)];
  yslice = [Y(i_best)];
  zslice = [Z(i_best)];
  slice(X,Y,Z,V,xslice,yslice,zslice)
  xlabel('\phi','FontSize', 14)
  ylabel('Tpot [C]','FontSize', 14)
  zlabel('gs [m]','FontSize', 14)

  phi=num2str(X(i_best));
  Tpot=num2str(Y(i_best));
  gs=num2str(Z(i_best));

  max_val=num2str(max(P_mod.Posterior_scaled(:)))
  title(['max(P)=',max_val,' at \phi=',phi,',Tpot=',Tpot,',gs=',gs], 'FontSize', 14)
  hcb=colorbar;
  ylabel(hcb,'P(T, \phi, gs | Vs)','FontSize', 14)
end
