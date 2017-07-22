% Solves 1D mcKenzie equations for compaction pressure (excess pressure)
%  Finite volume discretization - solves for pressure on cell centers using
%    permeabilities on cell edges. z = 0 and z = zmax are cell edges, so
%    ghost cells are set by calculating necessary ghost to fix P or dPdx at
%    boundary. Assumes uniform mesh spacing. 
%
%  See pgs 75-77 in lab notebook 3 (2011-06-30 to ...) for nondim,
%  discretization, BC info. 
% 
%  Calling Procedure
%   [P A src] = P_solve_gen(phi,phis,phiref,dz,np,nxi,n1,n2,BCv,solver)
%  Input
%   phi    nondimensional melt fraction array. phi(1) and phi(n+2) are 
%          ghost cells. 
%   phis   staggered melt fraction. phi(1) is at z=0, phi(n+1) is at z=zmax
%   phiref reference melt fraction
%   dz     cell spacing 
%   np     permeability exponent. k ~ phi^np
%   nxi    compaction viscosity exponent. xi ~ phi^nxi
%   n1     switch for 1-phi approximation.
%            n1=1 solves for 1-phi, n1=0 solves for 1-phi=1
%   n2     switch for 1+4/3phi approximation
%            n2=1 solves for 1+(4/3)(phi^-nxi) 
%            n2=1 solves for 1+(4/3)(phi^-nxi)=1
%   gz  switch for gravity, glfag = 0 for no gravity
%   BCv    Two element array with boundary values BCv(1) is pressure value
%          at z=0, BCv(2) is pressure value at z=zmax. 
%   solver flag for solver to use
%          1 direct on full matrix with A\b
%          2 direct on full matrix with inv(A)*b
%          3 direct on sparse matrix with A\b
%          4 direct on sparse matrix with inv(A)*b
% Output
%   P      excess (compaction) pressure. P(1), P(n+2) are ghost cells.
%   A      the matrix constructed for inversion
%   src    right hand side 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,A,src] = P_solve_gen(phi,phis,perms,mus,SurfCoff,dz,material,BCv,BCType)

%   load material settings
    al = material.alpha;
    nxi = material.nxi; 
    n1 = material.n1; 
    gz = material.gz;
    muf = material.mufo;
    muso = material.muso; 
    drhog = material.g*(material.rhos - material.rhof); 

%  mesh infos
   n = numel(phi); % 2 cells are ghosts          
   cdz = 1/(dz*dz);
     
%  material properties      
   k = perms; % permeability (cell edges)   
   mus = muso*exp(-al*phi); % shear viscosity (cell centers)  
%    mus(mus>1e22)=1e22; 
   xi = mus.*(1+4/3*phi.^nxi); % compaction viscosity (cell centers)
%    xi(xi>1e24)=1e24; 
%    xi(xi>1e26)=1e26;
   Fc = (1-phi).^n1; % cell centered solid fraction
   Fs = (1-phis).^n1; % staggered solid fraction
   kmuf = k./muf; % k / muf
   Buoy = kmuf .* drhog ; % Buoyancy driving force
   
% build generic entries 
%  mask 
   m1 = 2:n-1;
%  generic right hand side entries     
   src = zeros(n,1); 
   src(m1) = -gz*(Buoy(m1).*Fs(m1)-Buoy(m1-1).*Fs(m1-1))/dz;  
   Surf = -(SurfCoff(m1).*phi(m1+1)+SurfCoff(m1-1).*phi(m1-1) ...
           -phi(m1).*(SurfCoff(m1) + SurfCoff(m1-1)))/dz/dz; 
       
   src(m1) = src(m1) + Surf;    

%  generic diagonal and off-diagonal entries         
   dia = zeros(n,1); 
   diap = zeros(n,1);
   diam = zeros(n,1);
   dia(m1,1) = -1./xi(m1)-(kmuf(m1-1)+kmuf(m1)).*Fc(m1)*cdz;
   diap(m1+1,1) = kmuf(m1).*Fc(m1+1)*cdz;
   diam(m1-1,1) = kmuf(m1-1).*Fc(m1-1)*cdz;
   % NOTE: using spdiags later, so the off diagonal's storage, diap,diam
   % is offset because some of the buffering entries get chopped when
   % building the sparse matrix.

%  gotta rescale. 
   maxd = max(abs(dia(dia~=0))); 
   dia = dia./maxd; 
   diap = diap./maxd; 
   diam = diam./maxd; 
   src = src./maxd; 
   
% fill ghost cells
%  z = 0        
   if BCType(1) == 1 % dirichlet condition on P             
       dia(1) = 1;
       diap(2) = 1;
       src(1) = 2*BCv(1);
   elseif BCType(1) == 2  % neumann flux condition on P
       dia(1) = -1;
       diap(2) = 1;
       src(1) = BCv(1)*dz;
   end
%  z = zmax
   if BCType(2) == 1 % dirichlet condition on P                          
       dia(n) = 1;
       diam(n-1) = 1;
       src(n) = 2*BCv(2);
   elseif BCType(2) == 2 % neumann flux condition on P
       dia(n)= 1;
       diam(n-1)= -1;
       src(n) = BCv(2)*dz;
   end
   
% Assemble and solve     
  A = spdiags([diam dia diap],[-1 0 1],n,n);
%   disp(condest(A))
  P = A\src;
%   [L,U,P1]=lu(A); 
%   y = L\src; 
%   P = U\y;
%     P = pcg(A,src);
%   P = inv(A)*src; 
  
end
