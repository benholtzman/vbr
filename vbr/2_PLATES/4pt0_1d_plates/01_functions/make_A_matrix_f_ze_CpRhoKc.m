function [A] = make_A_matrix_f_ze_CpRhoKc(nu,Rho,Cp,K)
% Makes the A matrix operator to calculate gradients in T
% with spatially variable density, specific heat and conductivity
%
%
%    INPUT:
%      nu=dt/dx^2,
%      rho node-centered vector of density
%      cp node-centered vector of specific heat
%      K is a node-centered vector of conductivity
%
%    OUTPUT:
%      A(K,stencil for T)
%

%   make K a column vector and get matrix size
K=K(:); 
N=length(K);

%  calculate the conductivity at the half points - simple averaging
K_half=(K(1:N-1)+K(2:N))/2.;

%  Allocate space for sparse matrices 
A = spalloc(N,N,3*N); % discrete diffusion operator matrix
%  I = speye(N); % sparse identity matrix


%  Calculate interior rows for Diffusion Operator (Ku_x)_x
%  Loop for clarity of what's happening - implemented in efficient way
%  below using spdiags
% for i=2:N-1
%     j=i-1:i+1;
%     A(i,j)=nu*[ K_half(i-1)  -(K_half(i-1)+K_half(i))  K_half(i)]; 
% end

inv_RhoCp = 1./(Rho.*Cp);
K_half_m1=inv_RhoCp(2:N-1).*K_half(1:N-2); % minus 1
K_half_p1=inv_RhoCp(2:N-1).*K_half(2:N-1); % plus 1
K_half_c = -inv_RhoCp(2:N-1).*(K_half(1:N-2) + K_half(2:N-1)); %  center
 

A(2:N-1,:) = nu*spdiags([K_half_m1 K_half_c K_half_p1],0:2,N-2,N);

%A(2:N-1,:) = nu*spdiags([K_half(1:N-2) -(K_half(1:N-2) + K_half(2:N-1)) K_half(2:N-1)],0:2,N-2,N);


% BCs are set via ghosts outside of here... THIS IS ONLY OK FOR EXPLICIT FORWARD STEP. NEED
% TO MODIFY THE EQS IF YOU WANT TO INVERT A.
 A(1,:) = 0; 
 A(N,:) = 0; 

% Boundary Conditions (Dirichlet) - Constant Temperature
%A(1,1) = -2*K_half(1) ;
%A(1,2) = 2*K_half(1) ;
%A(N,N) = -2*K_half(end) ;
%A(N,N-1) = 2*K_half(end) ;

% % set left and right crank Nicolson matrices
% AL = I -  .5*A;
% AR = I + .5*A;

% %   set dirichlet boundary conditions for point i and N
% AL(1,1)=1; AR(1,1)=1;
% AL(N,N)=1; AR(N,N)=1;

end
