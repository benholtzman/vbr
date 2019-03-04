function [resids,best]=jointResiduals(var1name,pred_1,obs_1,var2name,pred_2,obs_2)
  % calculates joint residual for two variables
  % input:
  %   var1name   fieldname for final structure for 1st variable
  %   pred_1     prediction(s), matrix of any size for 1st variable
  %   obs_1      observed value, scalar for 1st variable
  %   var2name   fieldname for final structure for 2nd variable
  %   pred_2     prediction(s), matrix of any size for 2nd variable
  %   obs_2      observed value, scalar for 2nd variable
  %
  % See manual, Menke book Ch 11
  %
  % pred_1,pred_2 must be the same size
  [chi_sq_1,P1]=singleParamResids(obs_1,pred_1);
  [chi_sq_2,P2]=singleParamResids(obs_2,pred_2);
  P_Joint=(P1.*P2).^2;
  resids.(var1name)=chi_sq_1;
  resids.(var2name)=chi_sq_2;
  resids.P_Joint=P_Joint;

  best=struct();
  [~, i_best] = max(P_Joint(:));
  [var1_i, var2_i] = ind2sub(size(pred_1), i_best);
  best.('var1_i')=var1_i;
  best.('var2_i')=var2_i;

end

function [chi_sq,P1]=singleParamResids(obs,pred)
% single parameter residual components
  chi_sq=(pred - obs).^2./obs; % chi-squared
  Res_N = chi_sq ./ max(chi_sq(:));
  P1=(2*pi*Res_N).^-0.5 .* exp(-0.5*Res_N);
end
