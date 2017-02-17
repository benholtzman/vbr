%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pargasite_storage
%
% loads grabit points from Green et al., "Water and its influence on the 
% lithosphere-asthenosphere boundary," Nature 2010 fig 1b then interpolates
% to calculate the storage capacity of water in pargasite on the
% compuational grid. Then calculates an effective partition coefficient. 
%
% Input
%   z            the computational grid [km]
%   Cbulk        current bulk water concentration [wt %]
%   parg_file    points to the storage capacity data file
%   settings     settings structure from initialization 
%
% Output
%   Strg         storagce capcity [wt %]
%   Kd           effective partition coefficient 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Strg,Kd]=pargasite_storage(z,Cbulk,settings,Strg,loadparg)

if strcmp(loadparg,'load')
% get storage capacity on the computational grid
  [Strg,Coarse,PargStab]=interp_parg(z(2:end-1),settings.parg_file); 
  Strg=[Strg(1); Strg(1:end); Strg(end)]; 
end

Strg = Strg * strcmp(settings.Flags.parg,'yes');

% set initial partition coefficient vs depth -- will recalculate this so
% that Kd returns to H2O value when storage capacity is exceeded. 
  Kd = settings.kd_H2O;%*(Cbulk>=Strg) + settings.kd_parg*(Cbulk<Strg); 
%   disp(PargStab)
end

function [Fine,Pargasite,BreakdownZ]=interp_parg(z_fine,parg_file)
load(parg_file)

BreakdownZ=Pargasite(end,2); % the breakdown depth of pargasite
Z_Coarse=Pargasite(:,2); % coarse grid
Coarse_Val=Pargasite(:,1); % storage capacity
Fine=interp1(Z_Coarse,Coarse_Val,z_fine); 
Fine(z_fine>BreakdownZ)=0; 
end
