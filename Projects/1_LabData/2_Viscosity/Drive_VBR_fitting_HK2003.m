  close all; clear
  
%% add VBR paths 
   Work.VBR_version = 'VBR_v0p95';   
   addpath(['../../4_VBR/'  Work.VBR_version ],...
           ['../../4_VBR/'  Work.VBR_version '/functions'],...
           ['../../4_VBR/'  Work.VBR_version '/params'])
 

   
%%%%%%%%%%%%%%%%%%%%%%%    
%% Test the dry case %%
%%%%%%%%%%%%%%%%%%%%%%%

%% the thing to vary
   stress=logspace(1,3,55); 

%% construct state variable fields          
   VBR.in.SV.T_K = 1250 + 273;
   VBR.in.SV.P_GPa = 300 * 1e6/1e9;
   VBR.in.SV.sig_MPa = stress;
   VBR.in.SV.dg_um = 15;
   VBR.in.SV.phi = 1e-5;
   VBR.in.SV.Ch2o = 0;
   
%% write method list (these are the things to calculate)
   VBR.in.viscous.methods_list={'HK2003'};
   VBR.in.viscous.HK2003 = Params_Viscous('HK2003');
   VBR.in.viscous.HK2003 = rmfield(VBR.in.viscous.HK2003,'gbs');
%    VBR.in.viscous.HK2003.gbs.x_phi_c_gt1250 = 1;
%    VBR.in.viscous.HK2003.gbs.x_phi_c_lt1250 = 1;
%    VBR.in.viscous.HK2003.diff.x_phi_c = 1;
   
%% calculate viscosity   
   [VBR] = VBR_spine(VBR);  

%% load the grabit points
   load ./data/HK2003_fig2a_data_grab
   
%% plot   
   figure('color',[1 1 1])   
   loglog(sigma_Disl,sr_Disl,'.k','displayname','data grab: disl component')
   hold on 
   loglog(sigma_MK2000,sr_MK2000,'ok','displayname','data grab')
   
   loglog(stress,VBR.out.viscous.HK2003.diff.sr,'--k','displayname','diff')   
   loglog(stress,VBR.out.viscous.HK2003.disl.sr,'k','displayname','disl')
   if isfield(VBR.out.viscous.HK2003,'gbs')
       loglog(stress,VBR.out.viscous.HK2003.gbs.sr,'--b','displayname','gbs')
   end
   loglog(stress,VBR.out.viscous.HK2003.sr_total,'k','linewidth',2,...
          'displayname','composite')   
   xlabel('\sigma [MPa]')
   ylabel('\epsilon [s^{-1}]')
   title('Dry: see figure 2A of Hirth and Kohlstedt 2003')
   legend('location','northwest')
   
%%%%%%%%%%%%%%%%%%%%%%%    
%% Test the wet case %%
%%%%%%%%%%%%%%%%%%%%%%%

%% the thing to vary 
   Ch2o = logspace(1,4,100); 
   
%% construct state variable fields          
   VBR.in.SV.T_K = 1250 + 273;
   VBR.in.SV.P_GPa = 300 * 1e6/1e9;
   VBR.in.SV.sig_MPa = 150;
   VBR.in.SV.dg_um = 15;
   VBR.in.SV.phi = 0.0011;
   VBR.in.SV.Ch2o = Ch2o;
   
%% calculate viscosity   
   [VBR] = VBR_spine(VBR);  

%% load the grabit points
   load ./data/HK2003_fig5b_data_grab
   
%% plot   
   figure('color',[1 1 1])
   loglog(Ch2o_HK,sr_HK,'or','displayname','HK2003,fig5B')
   hold on
   loglog(Ch2o,VBR.out.viscous.HK2003.sr_total,'k','linewidth',2,...
          'displayname','composite')   
   loglog(Ch2o,VBR.out.viscous.HK2003.diff.sr,'--k','displayname','diff')
   loglog(Ch2o,VBR.out.viscous.HK2003.disl.sr,'k','displayname','disl')
   if isfield(VBR.out.viscous.HK2003,'gbs')
       loglog(Ch2o,VBR.out.viscous.HK2003.gbs.sr,'--b','displayname','gbs')
   end
   xlabel('C_{H_2O} [PPM]')
   ylabel('\epsilon [s^{-1}]')
   title('Wet: see figure 5B of Hirth and Kohlstedt 2003')
   legend('location','northwest')
      
   