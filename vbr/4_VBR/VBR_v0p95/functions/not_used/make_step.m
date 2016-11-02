function[step_vec] = make_step(N,N_z_steps,i_mid_step)
% make a step
% N = length of vector
% N_z_steps = length of 


x = linspace(-2,2,N_z_steps) ; 
step = (0.5.*(erf(x)+1))  ;



step_vec = (zeros(N,1)) ;   

top_trans = i_mid_step - (N_z_steps/2) ;
bot_trans = i_mid_step + (N_z_steps/2)-1  ;

count = 1 ;
for i = top_trans:(bot_trans-1)
    step_vec(i) = step(count) ;
    count = count+1 ;
end

for i = bot_trans:N
    step_vec(i) = 1 ;
end
