%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates weighted variance for a single parameter. 
%
% total variance of y, sig_y^2, is given by the sum of the
% partial-derivative weighted variance of each variable: 
% 
% sig_y^2 = sig_x1 ^ 2 * dy/dx1 + sig_x2 ^2 * dy/dx2 + ... 
% 
% where x1, x2 and ... are the independent parameters, each with a
% corresponding variance, sig_x1^2, sig_x2^2, ..., and partial derivative,
% dy/dx1, dy/dx2, ...
%
% this routine calculates just sig_x1 ^2 * dy/dx1. Can be called multiple
% times to generate suite of individual weighted variances to sum. 
%
% input 
%      y           values of function y at x
%      x           independent parameter array
%      stand_dev   the 1-sigma standard deviation of x 
%
% output
%      sig_x_square  the weighted-variance for x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function sig_x_square = single_param_variance(y,x,sig_x1)

% calculate the derivative of y with respect to x over full range using
% a central difference
  dydx=((y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2)));
  dydx_ave=mean(dydx);   
  % note: could add a check to see if there is a large variation in the
  % derivative over the range of interest. if there is, there might be more
  % error coming from a dependence of dydx on x. 
  
% calculate weighted-variance  
sig_x_square = sig_x1.^2 * (dydx_ave)^2; 

end