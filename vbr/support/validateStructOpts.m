function Options=validateStructOpts(func_name,func_varargin,Options,ValidOpts)
% general validator function for validating single-valued input options.
% overwrites valid fields in Options.
% validateStructOpts(func_name,func_vargin,Options,ValidOpts)
%  func_name = the name of the function calling this one
%  func_varargin = varargin from the function calling this one
%  Options = option structure with field name and default value
%  ValidOpts = structure with valid options for each field name


% handle the input arguments
  optionNames=fieldnames(Options);
  for option_pair = reshape(func_varargin,2,[])
   optName = lower(option_pair{1});
   optVal=option_pair{2};
   % check if this option is in list of valid options
   if any(strcmp(optName,optionNames))
      % check if the value of this option is in list of valid option values
      safeVals=ValidOpts.(optName);
      if any(strcmp(optVal,safeVals))
        Options.(optName) = optVal;
      else
        disp(['WARNING: in call to ',func_name])
        disp(['    ',optVal, ' is not valid for ',optName,' using default.'])
        disp(['    possible values are: ',strjoin(safeVals,' ,')])
      end
   else
      msg=['Call to ',func_name,' failed. ',optName,...
          ' is not a valid option. Valid options: ',strjoin(optionNames,' , ')];
      error(msg)
   end
  end

end
