function TestResults = run_tests(test_type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TestResults = run_tests(test_type)
%
% runs all the test functions
%
% Parameters
% ----------
% test_type    optional string declaring the test_type to run. Options are:
%              'full_test'  runs all .m functions in this directory.
%
%
% Output
% ------
% TestResults   Structure with test results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~exist('test_type','var')
    test_type='full_test';
  end

  disp(['Running test_type: ',test_type]); disp('')

  % initialize VBR
  addpath('../..')
  vbr_init

  % run test functions
  TestResults=struct();
  failedCount=0;
  disp(['Starting ',test_type]); disp(' ')
  if strcmp(test_type,'full_test')
    % runs all test functions in this directory
    mfiles=dir('*.m');
    for ifile = 1:numel(mfiles)
      fname=mfiles(ifile).name;
      if ~strcmp('run_tests.m',fname)
        [fdir,funcname,ext]=fileparts(fname);
        try
          testResult=feval(funcname);
          disp('    test passed :D'); disp(' ')
        catch
          disp(['    ',funcname,' failed :('])
          disp(['    please run ',funcname,'() and debug.']); disp(' ')
          testResult=false;
          failedCount=failedCount+1;
        end
        TestResults.(funcname)=testResult;
      end
    end
  else
    disp(['test_type ',test_type,' does not exist.'])
  end

  % display the failed test functions
  disp('Testing complete.')
  disp(' ')
  if failedCount>0
    disp('Displaying failed test functions. Please run each one and debug:')
    fldz=fieldnames(TestResults);
    for ifi = 1:numel(fldz)
      fld=TestResults.(fldz{ifi});
      if fld==0
        disp(['    ',fldz{ifi}])
      end
    end
  else
    disp('all test functions ran successfully')
  end

end
