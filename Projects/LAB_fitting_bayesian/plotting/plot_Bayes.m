function plot_Bayes(posterior, sweep, obs_name, q_method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot_Bayes(posterior, sweep, obs_name)
%
% Plot the posterior probability distribution across the parameter ranges
% in sweep.
%
% Parameters:
% -----------
%      posterior           (size(sweep.Box)) matrix of posterior
%                          probability for each parameter combination
%
%      sweep               structure with the following fields
%               state_names     cell of the names of the varied parameters
%               [param name]    vector of the range of values that were
%                               calculated
%               Box             output of VBR calculation
%               (also other fields recording values relevant to the
%               calculation)
%
%       obs_name            string of the observation name, e.g. Vs, Qinv
%                           Only used to label the figure.
%       
%       q_method            string of the method to use for calculating
%                           the anelastic effects on seismic properties
%                           Only used to label the figure.
%
% Output:
% -------
%       Plot of the posterior probability distribution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  figure()
  q_method = strrep(q_method, '_', '\_');
  
  max_val = max(posterior(:));
  i_best = find(posterior(:) == max_val);
  fields = sweep.state_names;
  
  gs_ind = find(strcmp('gs', fields));
  if ~isempty(gs_ind)
      sweep.(fields{gs_ind}) = sweep.(fields{gs_ind}) ./ 1e6;
  end
      
  params = make_param_grid(fields, sweep);

  Y = params.(fields{1});
  X = params.(fields{2});
  fields_str = [fields{1}, ', ', fields{2}];

  if length(fields) == 3
      Z = params.(fields{3});
      fields_str = [fields_str ', ', fields{3}];
      xslice = X(i_best);
      yslice = Y(i_best);
      zslice = Z(i_best);
      slice(X, Y, Z, posterior, xslice, yslice, zslice)
      zlabel(fields{3},'FontSize', 14)
      
      title({['max(p) = ', num2str(max_val),' at ']; ...
          [fields{2}, ' = ', num2str(X(i_best)), ', ' ...
          fields{1}, ' = ', num2str(Y(i_best)), ', ', ...
          fields{3} ' = ', num2str(Z(i_best))]; ...
          ['using ', q_method]}, 'FontSize', 14)
      
  else
      imagesc(X(1,:), Y(:,1), posterior);
      title({['max(P)=', num2str(max_val),' at ']; ...
          [fields{2}, ' = ', num2str(X(i_best)), ', ' ...
          fields{1}, ' = ', num2str(Y(i_best))]; ...
          ['using ', q_method]}, 'FontSize', 14)
  end
  
  cblabel = ['P(', fields_str, ' | ', obs_name ')'];
  xlabel(fields{2},'FontSize', 14)
  ylabel(fields{1},'FontSize', 14)
  
  hcb = colorbar;
  ylabel(hcb,cblabel,'FontSize', 14)
end