function [M, E_G, E_list] =  my_dmrg (MPO, N_keep, N_sweep, MPS_choice, update_rule, MPO_length, physical_dim) 
     if MPS_choice == 0
          M_init = my_iter_diag(MPO_length, MPO, N_keep);
     elseif MPS_choice == 1
          M_init = random_mps(MPO_length, physical_dim, N_keep);
     end

     if update_rule == 1
          [M, E_G, E_list] = my_1site (M_init, MPO, N_keep, N_sweep);
     elseif update_rule == 2
          [M, E_G, E_list] = my_2site (M_init, MPO, N_keep, N_sweep);
     end
end

%MPS choice: 0 for iterative diagonalization and 1 for random MPS
%update rule: 1 for 1 site update and 2 for 2 site update




