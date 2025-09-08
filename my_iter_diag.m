function M_init = my_iter_diag(L, Hs, N_keep)
M_init = cell(1,L);

%initialize as dummy tensor
H_prev = 1;
A_prev = 1;

physical_dim = size(Hs{1}, 1);

for i = (1:L)
    left_dim = size(A_prev, 2);
    right_dim = left_dim * physical_dim;
    iden_now = eye(right_dim);
    iden_now = reshape(iden_now, [left_dim, right_dim, physical_dim]);
    H = contract(H_prev, 3, 2, iden_now, 3, 1);
    H = contract(H, 4, [2, 4], Hs{i}, 4, [3, 2]);
    H = contract(H, 4, [1, 3], iden_now, 3, [1, 3], [3, 1, 2]);
    %Now diagonalize the Physical Hamiltonian of the system
    H_phy = H(:,:,1); %physical hamiltonian matrix we want
    [V, D] = eig((H_phy + H_phy')/2); %수치적 오류로 발생한 허수부 제거
    [D, idx] = sort(diag(D), 'ascend');
    %남길 차원 수 계산
    if i < L
        num_reserve = min([numel(D);N_keep]);
    else
        num_reserve = 1;
    end
    %필요한 고유 벡터만 남김
    V = V(:, idx(1:num_reserve));
    A_now = contract(iden_now, 3, 2, V, 2, 1, [1, 3, 2]);
    M_init{i} = A_now;
    A_prev = A_now;

    % basis change to V
    H_prev = contract(H, 3, 2, V, 2, 1);
    H_prev = contract(H_prev, 3, 1, conj(V), 2, 1, [3, 2, 1]);

end

end
