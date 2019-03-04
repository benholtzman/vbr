%%  ==================================================
%%  FUNCTIONS (now go last !)
%%  ==================================================


function [i_T_d1, i_g_d2, i_P_d3] = find_index_f(VBR,state)
  disp('data conditions:')
  state ;
  T_C = state.T_C
  T_C_vec = VBR.in.SV_vectors.T_K_vec_dim1 - 273 ;
  i_T_d1 = find(T_C_vec > T_C ,1)-1 ;
  T_C_vec(i_T_d1)

  g_um = state.dg_0
  gs_um_vec = VBR.in.SV_vectors.gs_um_vec_dim2 ;
  i_g_d2 = find(gs_um_vec > g_um ,1)-1 ;
  gs_um_vec(i_g_d2)

  P_GPa = state.P_GPa
  P_GPa_vec = VBR.in.SV_vectors.P_GPa_vec_dim3 ;
  i_P_d3 = find(P_GPa_vec > P_GPa ,1)-1 ;
  P_GPa_vec(i_P_d3)
end
