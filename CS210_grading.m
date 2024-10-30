clear
A=[1 1 0;1 2 2;0 2 2];
[L,U,P] = lu_pivoting(A);
L
U
P

function [L, U, P] = lu_pivoting(A)
    % LU Decomposition with Partial Pivoting
    n = size(A, 1);
    L = eye(n); % Initialize L as identity matrix
    U = zeros(size(A)); % Initialize U as zero matrix
    P = eye(n); % Initialize P as identity matrix
    A2 = A;

    for k = 1:n
        % Partial Pivoting
        [~, p] = max(abs(A2(k:n, k)));
        p = p + k - 1;
        if p ~= k
            % Swap rows in A2
            A2([k, p], :) = A2([p, k], :);
            % Swap rows in P
            P([k, p], :) = P([p, k], :);
            % Swap entries in L (columns 1 to k-1)
            % if k > 1
            %     L([k, p], 1:k-1) = L([p, k], 1:k-1);
            % end
        end
        if A2(k, k) == 0
            error('Zero pivot encountered');
        end
        % Compute U(k, k:n)
        U(k, k:n) = A2(k, k:n);
        % Compute L(k+1:n, k)
        L(k+1:n, k) = A2(k+1:n, k) / A2(k, k);
        % Update A2
        for i = k+1:n
            A2(i, k+1:n) = A2(i, k+1:n) - L(i, k) * U(k, k+1:n);
        end
    end
end
