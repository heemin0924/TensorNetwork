Heemin Kim, 2021-12752, Department of Physics and Astronomy 

- For Problem 1, I implemented my_dmrg.m with four configuration options.

When MPS_choice is set to 0, the algorithm uses an initial MPS obtained from iterative diagonalization.
When MPS_choice is set to 1, it uses a randomly generated MPS.
The update rule is determined by the update_rule input: 1 corresponds to the single-site update, while 2 triggers the two-site update scheme.

The code includes several custom functions: my_iter_diag.m, random_mps.m, my_1site.m, and my_2site.m, along with my_lanczos_1site.m and my_lanczos_2site.m used internally.
Additionally, I utilized the contract function and the svdTr function in the implementation.
The overall procedure is described in detail with pseudocode in the attached PDF, which corresponds to Problem 2.

- For Problem 3, I made tb_chain.mlx

I designed the tb_chain_MPO function to construct the MPO representation of a tight-binding chain. For each sub-problem, I input a different list of hopping coefficients. I then examined how the error between the ground-state energy obtained from my_dmrg and the exact result from nonIntTB evolves during the sweeps and as N_keep varies. This analysis was repeated for every combination of initial MPS and update method.

- For Problem 4, I made hubbard.mlx

I designed the Hubbard_MPO function to construct the MPO representation of a 1D Hubbard model. I examined how the ground-state energy after the sweeps behaves for various values of N_keep, and estimated the energy in the limit as N_keep approaches infinity by stopping when the ground-state energy showed negligible change. Then, I measured the double occupancy along the chain by computing it site by site, moving the orthogonality center through the chain using SVD. This procedure was repeated for the case of U = 10 (the first case used U = 2).



