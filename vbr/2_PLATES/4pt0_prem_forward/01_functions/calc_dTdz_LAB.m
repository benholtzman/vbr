%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the temperature gradient that enforces the given heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dTdz = calc_dTdz_LAB(Vark,settings)

%     T = Vark.T(end-1); % temperature [C]
    k = (Vark.Kc(end)+Vark.Kc(end-1))/2; % thermal conductivity [W/m/K]
    QLAB = settings.Q_LAB; % [W/m2]
    
    dTdz = - QLAB / k; 
    
    


end