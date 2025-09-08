function [M, E_G, E_list] = my_2site (M, H, N_keep, N_sweep)
    %variational energy list obtained over iterations
    L = numel(M); %site number
    E_list = zeros(L-1, 2*N_sweep);

    %hyperparameters
    nKrylov = 5;
    tol = 1e-8;

    %save effective Hamiltonian
    H_eff = cell(1,L+2);
    H_eff{1} = 1;
    H_eff{end} = 1;

    for i = 1:L
        X = contract(H_eff{i}, 3, 2, M{i}, 3, 1);
        Y = contract(X, 4, [2, 4], H{i}, 4, [3, 2]);
        H_eff{i+1} = contract(Y, 4, [1, 3], M{i}, 3, [1, 3], [3, 1, 2]);
    end

    for it = 1:N_sweep
        %right to left half-sweep
        for i = L-1:-1:1
            M_old = contract(M{i}, 3, 2, M{i+1}, 3, 1, [1, 3, 2, 4]);
            [M_new, E_list(i, 2*it-1)] = my_lanczos_2site(H_eff{i}, H{i}, H{i+1}, H_eff{i+3}, M_old, nKrylov, tol); %cite update and save variational energy via Lanczos method
            [M{i}, S, M{i+1}] = svdTr(M_new, 4, [1, 3], N_keep, 1e-8); %we can safely truncate small singular values in 2 site update
            %update left tensor
            M{i} = contract(M{i}, 3, 3, diag(S), 2, 1, [1,3,2]);

            %update right effective hamiltonian
            X = contract(H_eff{i+3}, 3, 1, M{i+1}, 3, 2);
            Y = contract(X, 4, [2, 4], H{i+1}, 4, [4, 1]);
            H_eff{i+2} = contract(Y, 4, [3, 1], M{i+1}, 3, [3, 2], [1, 3, 2]);
        end
        msg = sprintf('Sweep #%d of %d (right -> left) : Energy = %.7g', ...
              2*it -1, 2*N_sweep, E_list(1, 2*it -1));
        disptime(msg);
        %left to right half-sweep
        for i = 1:L-1
            M_old = contract(M{i}, 3, 2, M{i+1}, 3, 1, [1, 3, 2, 4]);
            [M_new, E_list(i, 2*it)] = my_lanczos_2site(H_eff{i}, H{i}, H{i+1}, H_eff{i+3}, M_old, nKrylov, tol);
            [M{i}, S, M{i+1}] = svdTr(M_new, 4, [1, 3], N_keep, 1e-8);
            M{i} = permute(M{i}, [1, 3, 2]);
            %update right tensor
            M{i+1} = contract(diag(S), 2, 2, M{i+1}, 3, 1);

            %update effective Hamiltonian
            X = contract(H_eff{i}, 3, 1, M{i}, 3, 1);
            Y = contract(X, 4, [2, 4], H{i}, 4, [3, 1]);
            H_eff{i+1} = contract(Y, 4, [1, 3], M{i}, 3, [1, 3], [1,3,2]);
        end
        msg = sprintf('Sweep #%d of %d (left -> right) : Energy = %.7g', ...
              2*it, 2*N_sweep, E_list(L-1, 2*it));
        disptime(msg);
    end

    E_G = E_list(L-1, 2*N_sweep);

end




