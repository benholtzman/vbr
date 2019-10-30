function [resids,best]=jointResiduals(var1name, pred_1, obs_1, sigma_1,...
    var2name, pred_2, obs_2 ,sigma_2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates joint residual for two variables
% input:
%   var1name   fieldname for final structure for 1st variable
%   pred_1     prediction(s), matrix of any size for 1st variable
%   obs_1      observed value, scalar for 1st variable
%   sigma_1    uncertainty, scalar for 1st variable
%   var2name   fieldname for final structure for 2nd variable
%   pred_2     prediction(s), matrix of any size for 2nd variable
%   obs_2      observed value, scalar for 2nd variable
%   sigma_2    uncertainty, scalar for 2nd variable
%
% See manual, Menke book Ch 11
%
% pred_1,pred_2 must be the same size
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get chi-squared and individual P for both
  [chi_sq_1,P1]=singleParamResids(obs_1, pred_1, sigma_1);
  [chi_sq_2,P2]=singleParamResids(obs_2, pred_2, sigma_2);

% calculate joint P
  P_Joint=(P1.*P2).^2;

% store in resid structure
  resids.([var1name,'_chisq'])=chi_sq_1;
  resids.([var2name,'_chisq'])=chi_sq_2;
  resids.(['P_',var1name])=P1;
  resids.(['P_',var2name])=P2;
  resids.P_Joint=P_Joint;

% store best i,j in best structure 
  best=struct();
  [~, i_best] = max(P_Joint(:));
  [var1_i, var2_i] = ind2sub(size(pred_1), i_best);
  best.('var1_i')=var1_i;
  best.('var2_i')=var2_i;

end

function [chi_sq,P1]=singleParamResids(obs, pred, sigma)
% single parameter residual components
  chi_sq=(pred - obs).^2./sigma.^2; % chi-squared
  P1=(2*pi*sigma.^2).^-0.5 .* exp(-0.5*chi_sq);
end
