function[V_rlx] = V_rlx_calc(Vu_vec,Q_mat)

ra = size(Vu_vec) ;
nx1 = ra(1) ;
dims = length(ra);
if ra(2)==1
    dims=1    ;
elseif ra(2)>1
    nx2 = ra(2) ;
    if dims==3
        nx3 = ra(3) ;
    end
end

Qsize = size(Q_mat) ;
nf = Qsize(1) ; % (number of frequency values)

% MAKE matrixes with added dimension for frequency in (i). 
for i = 1:nf
    if dims==3
        Vu_mat(i,:,:,:) = Vu_vec(:,:,:) ;
    elseif dims==2
        Vu_mat(i,:,:) = Vu_vec(:,:) ;
    elseif dims==1
        Vu_mat(i,:) = Vu_vec(:) ; 
    end
end


% WHERE IS THIS FROM ???
alf = 0.2 ; 
V_rlx = Vu_mat.*(1-(2./(Q_mat.*tan(alf*pi/2)))) ;


%K_vec = RA.K_vrh_90.*1e9 ; % 
%rho_vec = RA.rho_eos ; 

%RA.Andrade_f(i).Vs_r = Vs_r ;
    %RA.Andrade_f(i).Vp_r = Vp_r ;
    
%        Vs_r = sqrt(G_r./(rho_vec))  ;
    %Vp_r = sqrt((K_vec+(4/3).*G_r)./(rho_vec)) ;