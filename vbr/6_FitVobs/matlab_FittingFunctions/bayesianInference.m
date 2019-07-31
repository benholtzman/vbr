function [probs] = bayesianInference(States,states_fields,Obs,Residuals, ...
    obs_field, ifnormal)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculate posterior probability distribution using Bayesian inference
  % assuming normal distribution of data and prior-model parameters.
  %
  % Bayes' theorem:  P(A|D) = ( P(D|A) P(A) ) / P(D)
  %   P(D|A) : likelihood (probability of measurements D occuring given A)
  %   P(A) : prior-model
  %   P(D) : probability of data
  %   P(A|D) : posterior probability of A given D.
  %
  % For example, with Vs = f(phi, Tp, gs),
  %   A is [phi, Tp, gs];  B is observed Vs
  % P(T, phi, gs | Vs) = P( Vs | T, phi,gs) * P(T,phi,gs) / P(Vs)
  %
  % Input:
  %    States    :   structure with thermodynamic states
  %          .var1 : values of var1 at this state e.g., T
  %          .var1_mean : mean of var1
  %          .var1_std : standard devation of var1
  %          .var2 : e.g., phi
  %          .var2_mean : mean of var2
  %          .var2_std : standard devation of var2
  %
  %          all vari variables must matrices of the same size:
  %                   States.Tpot(n_Tpot,n_phi,n_gs)
  %                   States.phi(n_Tpot,n_phi,n_gs)
  %                   States.gs(n_Tpot,n_phi,n_gs)
  %
  %    states_fields : cell array with base fieldnames in States, e.g., {'Tpot','phi'}
  %                   each fieldname must exist in States.
  %
  %    Obs       :    structure of observations
  %       .obs_field       :  observations e.g., Obs.meanVs=V1
  %       .obs_field_std   :  std of obs
  %       .obs_field_mean  :  mean of obs
  %
  %    Residuals : array-structure of chi-square residuals, same size as States.
  %       .obs_field
  %
  %    obs_field  : the Obs and Residuals field to use for inference e.g., 'meanVs'
  %
  % Ouput
  %    probs      :  structure with probability distributions
  %         .Prior_mod  : multivariate prior model PDF, P(A)=P(v1)*P(v2)...
  %         .P_Obs_given_mod : likelihood PDF, P(D|A)
  %         .P_Obs  : data PDF, P(D)
  %         .Posterior : unscaled posterior distrubtion
  %         .Posterior_scaled : scaled posterior distrubtion
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % observed probability of data P(D), e.g., P(Vs). Scalar value.
  % This is the probability of the observed Vs being correct
  % assumes a normal distribution, mean = observed Vs, std = obs error
  % sigmaObs = mean(Obs.([obs_field,'_error']); % std = mean error
  % P_Obs = 1/sqrt(2*pi*sigmaObs^2); % assume mean = observed Vs so (x-mu) = 0
  sigmaObs =  Obs.([obs_field,'_std']); % standard deviation
  mu     = Obs.([obs_field,'_mean']); % mean value
  x      = Obs.([obs_field]); % measurement
  probs.P_Obs = normpdf(x,mu,sigmaObs); % probability of Data.
  disp(['P(D):'])
  disp(probs.P_Obs)

  % calcualte prior model for states: P(var1,var2,...)=P(var1)*P(var2)*...
  [probs.Prior_mod,sigmaPreds]=priorModelProbs(States,states_fields,ifnormal);

  % calculate likelihood P(D|A), e.g., P(Vs | T, phi, gs), using residual.
  % given residual, P(D|A) = 1/sqrt(2*pi*residual) * exp(-residual/2)
  % See manual, Menke book Ch 11
  % residual(k) here is a chi-squared residual. Given chi-square, PDF of data
  % with normal distribution is
  %   P = 1/sqrt(2*pi*sigma^2) * exp(-0.5 * chi-square)
  % where sigma = std of data, chi-square=sum( (x_obs - x_preds)^2/sigma^2 )
  % e.g., https://www-cdf.fnal.gov/physics/statistics/recommendations/modeling.html
  probs.P_Obs_given_mod = (2*pi*sigmaObs^2).^-0.5 .* exp(-0.5*Residuals.(obs_field))

  % calculate Posterior distribution
  probs.Posterior = (probs.P_Obs_given_mod .* probs.Prior_mod)./(probs.P_Obs)

  % scale Posterior distribution
  norm_for_residual = (2*pi*sigmaObs.^2).^-0.5 .* exp(-0.5*sigmaObs.^2); % ???
  probs.Posterior_scaled = probs.Posterior .* (2*pi*sigmaPreds./sigmaObs) / norm_for_residual

end

function [Prior_mod, sigmaPreds] = priorModelProbs(States,states_fields,ifnormal)
  % loops over the fields of States, calculates pdf for each field and total
  % prior model pdf:
  %   P(var1,var2,...)=P(var1)*P(var2)*...

  sigmaPreds=1;
  Prior_mod=1;
  for i_field = 1:numel(states_fields)
    this_field=states_fields{i_field};
    std_field=[this_field,'_std']; % e.g., Tpot_std
    mn_field=[this_field,'_mean']; % e.g., Tpot_mean
    if ifnormal
      sigma =  States.(std_field); % standard deviation
      mu     = States.(mn_field); % mean value
      x      = States.(this_field); % measurements
      P_var_i = normpdf(x,mu,sigma); % 1/sqrt(2*pi*sigma^2) * exp((-(x-mu)^2)/(2*sigma^2));
    else
      sigma = 1;
      min_val = min(States.(this_field)(:));
      max_val = max(States.(this_field)(:));
      x = States.(this_field); % measurements
      P_var_i = unifpdf(x, min_val, max_val); % uniform PDF over total range 
    end
    sigmaPreds=sigmaPreds.*sigma;
    Prior_mod=Prior_mod.*P_var_i; % propagate the probability
  end
end
