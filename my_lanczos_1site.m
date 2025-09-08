function [M_new, E_G] = my_lanczos_1site(H_left, H_center, H_right, M_old, nKrylov, tol)

v0 = M_old(:) / sqrt(abs(dot(M_old(:), M_old(:))));

V = zeros(length(v0), nKrylov);
alpha = zeros(nKrylov, 1);
beta = zeros(nKrylov - 1, 1);

V(:,1) = v0;

for i = 1:nKrylov
    v = reshape(V(:,i), size(M_old));
    X = contract(H_left, 3, 2, v, 3, 1);
    Y = contract(X, 4, [2, 4], H_center, 4, [3, 2]);
    Hv = contract(Y, 4, [2, 4], H_right, 3, [2, 3], [1, 3, 2]);
    w = Hv(:);
    alpha(i) = real(dot(V(:,i), w)); 
    if i < nKrylov 
        for it = 1:2 % repeat 2 times to reduce numerical noise
            x = zeros(size(w));
            for j = 1:i
                x = x + dot(V(:,j), w) * V(:,j);
            end
            w = w - x;
        end
        beta(i) = sqrt(abs(dot(w, w)));
        if beta(i) < tol
            V = V(:,1:i);
            alpha = alpha(1:i);
            beta = beta(1:i-1);
            break;
        end
        V(:,i+1) = w / beta(i);
    end
end

T = diag(alpha) + diag(beta,1) + diag(beta,-1);
[W, D] = eig((T + T')/2); %reduce numerical error
[Emin, index] = min(diag(D));
vmin = V(:,1:length(W(:,index))) * W(:,index);
M_new = reshape(vmin, size(M_old));

X = contract(H_left, 3, 2, M_new, 3, 1);
Y = contract(X, 4, [4, 2], H_center, 4, [2, 3]);
Z = contract(Y, 4, [2, 4], H_right, 3, [2, 3]);
E_G = real(contract(conj(M_new), 3, [1,2,3], Z, 3, [1,3,2]));

end
