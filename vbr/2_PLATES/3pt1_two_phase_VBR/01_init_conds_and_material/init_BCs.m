%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       init_BCs                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    %
% [BCs]=init_BCs(Info,settings)                                      %
%                                                                    %
% where BCs is a structure with a field for each variable that needs %
% a boundary condition:                                              %
%   BCs.type_i(1:2) sets the BC type for the top (1) and bottom (2), %
%                   where the type can be:                           %
%                     1 for dirichlet with specified value           %
%                     2 for neumann fixed value                      %
%                     3 for time-dependent neumann where the normal  %
%                           derivative is continous at the boundary  %
%   BCs.val_i(1:2)  sets the BC value for the top (1) and bottom (2).%
%                                                                    %
%   The variables that need boundary conditions:                     %
%             excess pressure    P                                   %
%             melt fraction      phi                                 %
%             temperature        T                                   %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [BCs]=init_BCs(Info,settings)

%%%%%%%%%%%%%%%%%%%%%%
%   excess pressure
%   A few notes:
%    The excess pressure should be set so that the desired flux S is 
%    obtained. 
%
%    Prescribed flux: 
%            S=0 = -k/mu * (dPdz + deltaRho * g)
%    For an impermeable boundary, S = 0, and so the gradient in excess
%    pressure must balance the buoyant pressure gradient: 
%      S=0 = -k/mu * (dPdz + deltaRho * g)  -->  dP/dz = - deltaRho * g
%    There are two types of "open" boundaries possible: 
%      compaction-free boundaries, where P = 0
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_P(2) = 0; 
     BCs.type_P(2) = 1;  
     BCs.val_P(1) = 0; 
     BCs.type_P(1) = 2;      
%      BCs.val_P(1) = -(settings.rhos-settings.rhof)*settings.g*settings.gz; 
     
%%%%%%%%%%%%%%%%%%%%%%     
%   melt fraction  
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_phi(1) = 0;   
     BCs.type_phi(1) = 2;   
     BCs.val_phi(2) = Info.init.phi(end);%settings.phimin;   
     BCs.type_phi(2) = 1;   
     
%%%%%%%%%%%%%%%%%%%%%%     
%   temperature      
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_T(2) = Info.init.T(end)+(settings.Tpot_excess-settings.Tpot); 
     BCs.val_T(1) = Info.init.T(1); 
     BCs.type_T(2) = 1;   
     BCs.type_T(1) = 1; 

%%%%%%%%%%%%%%%%%%%%%%     
%   water concentration
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_Cs(1) = 0;   
     BCs.type_Cs(1) = 2;   
     BCs.val_Cs(2) = settings.Cs0;   
     BCs.type_Cs(2) = 1;   
     
%%%%%%%%%%%%%%%%%%%%%%     
%   CO2 concentration
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_Cs_CO2(1) = 0;   
     BCs.type_Cs_CO2(1) = 2;   
     BCs.val_Cs_CO2(2) = settings.Cs0_CO2;   
     BCs.type_Cs_CO2(2) = 1;   

%%%%%%%%%%%%%%%%%%%%%%     
%   general material property
%%%%%%%%%%%%%%%%%%%%%%

     BCs.val_gen(1:2) = 0;   
     BCs.type_gen(1:2) = 2;   
       
end
