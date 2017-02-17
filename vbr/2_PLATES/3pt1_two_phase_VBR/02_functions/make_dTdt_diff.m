function [dTdt] = make_dTdt_diff(z,Rho,Cp,K,T)
% finite difference calculation for diffusion with variable coefficients

nz = numel(z);
dTdt = zeros(nz,1); 

for iz = 2:nz-1
%   cell edge locations 
    z_m = (z(iz) + z(iz-1))/2; 
    z_p = (z(iz) + z(iz+1))/2; 
    
%   conductivity on cell edges    
    k_m = (K(iz-1)+K(iz))/2;
    k_p = (K(iz+1)+K(iz))/2;
    
%   flux on cell edges
    f_m = k_m * (T(iz)-T(iz-1))/(z(iz)-z(iz-1)); 
    f_p = k_p * (T(iz+1)-T(iz))/(z(iz+1)-z(iz)); 
    
%   now we calculate!
    dTdt(iz) =  1/(Rho(iz)*Cp(iz)) * (f_p - f_m)/(z_p - z_m);     
end


end
