function[VBR]=Q_Andrade_f(VBR) ;
%% =============================================
% SORT what comes in: eta_ss and Gu are the viscosity and unrelaxed modulus
% arrays that depend on the state variables... for now the state variables
% do not need to enter this calculation directly. 
% this could change ! 

f_vec = VBR.ISV.f ;
w_vec = 2*pi.*f_vec ; 



Mu_in = VBR.Gu_vec ; % Pa unrelaxed shear modulus VECTOR
Ju_in = 1./Mu_in ; 
rho_in = VBR.ISV.rho ;
%Vsu_in = VBR.Vs_SnLG_Unrlx ; 

%% Andrade parameters, set in Params file
n = VBR.Andrade.n  ; % andrade exponent 

if VBR.Andrade.scaling_opt==1
    A = VBR.Andrade.A ;    
    eta_ss = VBR.Andrade.eta_ss ; 
    eta_ss_in = eta_ss.*ones(size(Mu_in)) ;
elseif VBR.Andrade.scaling_opt==2
    % MIGHT BE THAT WE USE TWO VISCOSITIES -
    % ONE FOR GB sliding and one for TOTAL !! check later ! 
    eta_ss_in = VBR.eta_diff ; % Pa.s
    %eta_ss = VBR.eta_diff ; % Pa.s
    %A = VBR.Andrade.Fudge.*(1/n) .* (eta_ss.^-n) .* (Gu./3).^(n-1) ; 
    %A = (1/n) .* (eta_ss_mat.^-n) .* (Mu_mat./3).^(n-1) ; 
end

%% Make the matrixes ===========================================
sz_ra = size(Mu_in) ;
nx1 = sz_ra(1) ;
dims = length(sz_ra);
if sz_ra(2)==1
    dims=1    ;
elseif sz_ra(2)>1
    nx2 = sz_ra(2) ;
    if dims==3
        nx3 = sz_ra(3) ;
    end
end
n_freq = length(f_vec) ;

%eta_ss = eta_ss.*ones(sz_ra) ; 

% MAKE matrixes with added dimension for frequency in (i). 
for i = 1:n_freq
    if dims==3
        eta_ss_mat(i,:,:,:) = eta_ss_in(:,:,:) ;
        Ju_mat(i,:,:,:) = Ju_in(:,:,:) ; 
        rho_mat(i,:,:,:) = rho_in(:,:,:) ;
        %Vs_mat(i,:,:,:) = Vsu_in(:,:,:) ; 
    elseif dims==2
        eta_ss_mat(i,:,:) = eta_ss_in(:,:) ; 
        Ju_mat(i,:,:) = Ju_in(:,:) ; 
        rho_mat(i,:,:) = rho_in(:,:) ;
        %Vs_mat(i,:,:) = Vsu_in(:,:) ; 
    elseif dims==1
        eta_ss_mat(i,:) = eta_ss_in(:) ; 
        Ju_mat(i,:) = Ju_in(:) ; 
        rho_mat(i,:) = rho_in(:) ;
        %Vs_mat(i,:) = Vsu_in(:) ; 
    end
end
Mu_mat = 1./Ju_mat ; 

% for andrade only
if dims==1
    w_mat = zeros(n_freq,nx1) ; 
    for x1 = 1:nx1
        w_mat(:,x1) = w_vec(:) ;
    end
elseif dims==2
    w_mat = zeros(n_freq,nx1,nx2) ; 
    for x1 = 1:nx1
        for x2 = 1:nx2
            w_mat(:,x1,x2) = w_vec(:) ;
        end
    end
elseif dims==3 ;
    w_mat = zeros(n_freq,nx1,nx2,nx3) ; 
    for x1 = 1:nx1
        for x2 = 1:nx2
            for x3 = nx3
                w_mat(:,x1,x2,x3) = w_vec(:) ;
            end
        end
    end
end
    


%% THE ANDRADE MODEL: ======================================

if VBR.Andrade.scaling_opt==2
    Mu_mat = 1./Ju_mat ;
    A = (1/n) .* (eta_ss_mat.^-n) .* (Mu_mat./3).^(n-1) ; 
end


J1 = Ju_mat + A.*gamma(1+n).*(w_mat.^-n)*cos(n*pi/2) ; 
J2 = A.*gamma(1+n).*(w_mat.^-n).*sin(n*pi/2) + 1./(eta_ss_mat.*w_mat) ; 
    
Qa = J1./J2 ; 
Ma = (J1.^2 + J2.^2).^(-1/2) ; % relaxed shear modulus ; 
% in the future calculate relaxed VELOCITY here too ! 

%% THE elastic-GBS relaxation peak: ======================================
% bump = VBR.Andrade.bump  ;
Te = VBR.Andrade.Te ; % = 0.1 ; % not sure what this is
Tgbs = VBR.Andrade.Tgbs ; % = 0.1 ;% sec
Delta = VBR.Andrade.Delta ;% = 0.43 ; % Relaxation strength

%J1_peak(i,:) = Ju*Delta./(1+omega.^2*Tgbs^2) ;
%J2_peak(i,:) = Ju*Delta*omega*Te./(1+omega.^2*Tgbs^2) ;

%J1_gbs = Ju_mat.*(1+(Delta)./(1+Tgbs.^2.*w_mat.^2)) ; 
J1_gbs = Ju_mat.*Delta./(1+Tgbs.^2.*w_mat.^2) ;
J2_gbs = Ju_mat.*Delta.*(w_mat.*Te)./(1+Tgbs.^2.*w_mat.^2) ;
Qgbs = J1_gbs./J2_gbs ; 
Mgbs = (J1_gbs.^2 + J2_gbs.^2).^(-1/2) ;

J1_comp = J1 + J1_gbs ;
J2_comp = J2 + J2_gbs ;
Qcomp = J1_comp./J2_comp ; 
Mcomp = (J1_comp.^2 + J2_comp.^2).^(-1/2) ;

Va = sqrt(Ma./rho_mat); 
Va_comp = sqrt(Mcomp./rho_mat) ; 

%% BOOKKEEPING ============================================

%VBR.Andrade.Tau_ve = Tau_ve ; 
VBR.Andrade.A = A ;

VBR.Andrade.J1 = J1 ;
VBR.Andrade.J2 = J2 ;
VBR.Andrade.Qa = Qa ;
VBR.Andrade.Ma = Ma ; 
VBR.Andrade.Va = Va ;

VBR.Andrade.J1_gbs = J1_gbs ;
VBR.Andrade.J2_gbs = J2_gbs ;
VBR.Andrade.Qgbs = Qgbs ;
VBR.Andrade.Mgbs = Mgbs ; 

VBR.Andrade.J1_comp = J1_comp ;
VBR.Andrade.J2_comp = J2_comp ;
VBR.Andrade.Qcomp = Qcomp ;
VBR.Andrade.Mcomp = Mcomp ; 
VBR.Andrade.Va_comp = Va_comp ;
    
