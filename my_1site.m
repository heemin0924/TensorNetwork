function [M, E_G, E_list] = my_1site (M, H, N_keep, N_sweep)
    %variational energy list obtained over iterations
    L = numel(M); %site number
    E_list = zeros(L, 2*N_sweep);

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
        for i = L:-1:1
            [M{i}, E_list(i, 2*it-1)] = my_lanczos_1site(H_eff{i}, H{i}, H_eff{i+2}, M{i}, nKrylov, tol); %cite update and save variational energy via Lanczos method
            [U, S, M{i}] = svdTr(M{i}, 3, 1, N_keep, 0);
            US = contract(U, 2, 2, diag(S), 2, 1);
            %update left tensor
            if i > 1      
                M{i-1} = contract(M{i-1}, 3, 2, US, 2, 1, [1, 3, 2]);
            else
                M{i} = contract(US, 2, 2, M{i}, 3, 1); %singular value is 1
            end

            %update right effective hamiltonian
            X = contract(H_eff{i+2}, 3, 1, M{i}, 3, 2);
            Y = contract(X, 4, [2, 4], H{i}, 4, [4, 1]);
            H_eff{i+1} = contract(Y, 4, [3, 1], M{i}, 3, [3, 2], [1, 3, 2]);
        end
        msg = sprintf('Sweep #%d of %d (right -> left) : Energy = %.7g', ...
              2*it -1, 2*N_sweep, E_list(1, 2*it -1));
        disptime(msg);
        %left to right half-sweep
        for i = 1:L
            [M{i}, E_list(i, 2*it)] = my_lanczos_1site(H_eff{i}, H{i}, H_eff{i+2}, M{i}, nKrylov, tol);
            [M{i}, S, V] = svdTr(M{i}, 3, [1, 3], N_keep, 0);
            M{i} = permute(M{i}, [1, 3, 2]);
            SV = contract(diag(S), 2, 2, V, 2, 1);
            if i < L
                M{i+1} = contract(SV, 2, 2, M{i+1}, 3, 1);
            else
                M{i} = contract(M{i}, 3, 2, SV, 2, 1,[1,3,2]);
            end

            %update effective Hamiltonian
            X = contract(H_eff{i}, 3, 1, M{i}, 3, 1);
            Y = contract(X, 4, [2, 4], H{i}, 4, [3, 1]);
            H_eff{i+1} = contract(Y, 4, [1, 3], M{i}, 3, [1, 3], [1,3,2]);
        end
        msg = sprintf('Sweep #%d of %d (left -> right) : Energy = %.7g', ...
              2*it, 2*N_sweep, E_list(L, 2*it));
        disptime(msg);
    end

    E_G = E_list(L, 2*N_sweep);

end




