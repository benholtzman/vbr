function [A,Fsrc] = make_A_matrix_QLAB(dz,Rho,Cp,K,QLAB,zLAB,zLABid,z_m)
% Makes the A matrix operator to calculate gradients in T
% with spatially variable density, specific heat and conductivity, with 
% an interior boundary condition where the heat flux at zLAB (node zLABid)
% is fixed by QLAB. 
%
%
%    INPUT:
%      dz
%      rho node-centered vector of density
%      cp node-centered vector of specific heat
%      K is a node-centered vector of conductivity
%      QLAB interior heat flux condition
%      zLAB location to apply heat flux
%      zLABid node corresponding to zLAB
%
%    OUTPUT:
%      A(K,stencil for T)
%
    nu = 1/dz/dz; 
%   make K a column vector and get matrix size
    K=K(:); 
    N=length(K);

%  calculate the conductivity at the half points - simple averaging
   K_half=(K(1:N-1)+K(2:N))/2.;

%  Allocate space for sparse matrices 
   A = spalloc(N,N,3*N); % discrete diffusion operator matrix

%  construct the operator for generic entries
   inv_RhoCp = 1./(Rho.*Cp);
   K_half_m1=inv_RhoCp(2:N-1).*K_half(1:N-2); % minus 1
   K_half_p1=inv_RhoCp(2:N-1).*K_half(2:N-1); % plus 1
   K_half_c = -inv_RhoCp(2:N-1).*(K_half(1:N-2) + K_half(2:N-1)); %  center
   
   A(2:N-1,:) = nu*spdiags([K_half_m1 K_half_c K_half_p1],0:2,N-2,N);

% interior flux BC
%   create the RHS 
    Fsrc=zeros(N,1);   
%   Zero the rows in the sparse matrix
    A(zLABid,:)=0; 
    A(zLABid+1,:)=0;
%     A(zLABid+2:N,:)=0; 
%   Modify appropriate entries 
%   for the point above the boundary
    dzLAB1 = zLAB - z_m(zLABid);
    dzLAB = dz/2 + dzLAB1; 
    A(zLABid,zLABid) = -inv_RhoCp(zLABid)*K_half(zLABid)/dzLAB/dz; 
    A(zLABid,zLABid-1) = inv_RhoCp(zLABid)*K_half(zLABid-1)/dzLAB/dz;
    Fsrc(zLABid) = inv_RhoCp(zLABid+1)*QLAB/dzLAB;
%   for the point below  
    dzLAB1 = z_m(zLABid+1) - zLAB ;
    dzLAB = dz/2 + dzLAB1; 
    A(zLABid+1,zLABid+1) = -inv_RhoCp(zLABid+1)*K_half(zLABid+1)/dzLAB/dz; 
    A(zLABid+1,zLABid+2) = inv_RhoCp(zLABid+1)*K_half(zLABid+2)/dzLAB/dz;
    Fsrc(zLABid+1) = -inv_RhoCp(zLABid+1)*QLAB/dzLAB;
    
% Global BCs are set via ghosts outside of here... 
% THIS IS ONLY OK FOR EXPLICIT FORWARD STEP. 
% NEED TO MODIFY THE EQS IF YOU WANT TO INVERT A.
  A(1,:) = 0; 
  A(N,:) = 0; 


end
