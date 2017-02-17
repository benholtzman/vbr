%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LABInfo = find_LAB(Vark,z,settings,LABInfo) 
% solidus zLAB
           zLAB = max(z);zLABid=numel(z); 
           izLAB = 3; LAB_success = 'failed'; 
           while izLAB < numel(Vark.T)
               if Vark.T(izLAB) > Vark.Tsol(izLAB)% && Vark.phi(izLAB) > 1e-5;%settings.phimin
                   zLABid = izLAB; % first one below boundary
                   izLAB = numel(Vark.T) * 2; 
                   zLAB = z(zLABid); 
                   LAB_success = 'succeeded'; 
               else
                   izLAB = izLAB + 1;
               end
           end  
% vicsous zLAB
%          ratio based on average eta below
             viscLABid='looking'; iz = 5; nz = numel(Vark.T); 
             while strcmp(viscLABid,'found')==0 && iz < nz-1;
                 zrange=[nz-iz:nz];
                 eta = Vark.eta(zrange);
                 zeta = z(zrange); 
                 eta_ave=cumtrapz(zeta,eta)./(max(zeta)-min(zeta)); 
                 eta_ave=eta_ave(end);
                 
                 
                 if eta(1)/eta_ave>=10; 
                     viscLABid='found'; 
                     zLABeta=z(nz-iz); 
                     zLABetaid=nz-iz;
                 end
                 iz = iz+1; 
             end
             if strcmp(viscLABid,'looking')==1
                 disp('viscous LAB not found! setting equal to solidus LAB')
                 zLABeta=zLAB; 
                 zLABetaid=zLABid;
             end
% melting onset
           
           izMO = numel(Vark.T)-2; 
           zMOid = izMO; zMO=z(zMOid);
           while izMO > 5
               if Vark.T(izMO) > Vark.Tsol(izMO)% && phi(izLAB) > settings.phimin
                   zMOid = izMO; % first one below boundary
                   zMO=z(zMOid); 
                   izMO = 1;
               else
                   izMO = izMO - 1;
               end
           end  
%    store it all	     
     LABInfo.zMOid=zMOid;
     LABInfo.zMO=zMO;
     LABInfo.zLABeta=zLABeta; 
     LABInfo.zLABetaid=zLABetaid; 
     LABInfo.zLAB=zLAB;
     LABInfo.zLABid=zLABid;
     LABInfo.phiLAB=Vark.phi(zLABid);
end
