function M_init = random_mps(L, d, N_keep)

M_init = cell(1,L);

left_dim = 1;
physical_dim = d;
right_dim = min(d, N_keep);

for i = 1:L
    if i == L
        right_dim = 1;
    end
    
    [A, ~] = qr(randn(left_dim * physical_dim, right_dim), 0);
    A = reshape(A, [left_dim, right_dim, physical_dim]);
    M_init{i} = A;
    left_dim = right_dim;
    if right_dim * physical_dim < N_keep
        right_dim = right_dim * physical_dim;
    else
        right_dim = N_keep;
    end
end
end




