function init_bilevel_toolbox()
%INIT_UNLOCBOX Initialize the toolbox
%   Usage: init_bilevel_toolbox()

  % Adding dependencies
  global GLOBAL_path;
  GLOBAL_path = fileparts(mfilename('fullpath'));

  addpath(genpath(GLOBAL_path));

  % Load version number
  bp = [GLOBAL_path,filesep];
  [FID,MSG] = fopen([bp,'bilevel_toolbox_version'],'r');
  if FID == -1
    error(MSG)
  else
    bilevel_toolbox_version = fgetl(FID);
    fclose(FID);
  end
