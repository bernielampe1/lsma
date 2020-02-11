function [A, S] = VNMF(HIM, c)
    k = 0;      % iter count
    maxiter = 10; % termination
    alpha = 0.001;  % A step size
    beta = 0.001;   % S step size
    lambda = 100; % volume parameter
    delta = 20; % ASC parameter
    
    % derivative parameters
    tao = lambda/factorial(c-1);
    C = [ones(1, c); zeros(c-1, c)];
    B = [zeros(1, c-1); eye(c-1)];

    % prep data
    X = linearize_bands(HIM)';
    [~, N] = size(X);
    X_tilde = [X; delta * ones(1, N)];
    
    % PCA of X
    mu = mean(X, 2);
    K = X - repmat(mu, 1, N);
    [U, ~] = eig(K*K');
    U = fliplr(U);
    U = U(:,1:(c-1));
    
    % abundance fractions init
    S = zeros(c, N);
    
    % init A endmember matrix with data
    A = X(:, randi([1, N], 1, c));
    
    % loop
    while k < maxiter
        % learning rule for A
        Z = C + B * U' * (A - mu * ones(1, c));
        delA = (A * S - X) * S' + 1000 * det(Z).^2 * U * B' * inv(Z)';

        A = A - alpha * delA;
%        A(A<0) = 0;

        % augment A for ASC and learning rule for S
        %A_tilde = [A; delta * ones(1, c)];
        %S = S - beta * (A_tilde' * (A_tilde * S - X_tilde));
        S = S - beta * (A' * (A *S - X));
 %       S(S<0) = 0;
        
        % scale
 %       S = normc(S);
        
        % iter increment
        k = k + 1;
    end
end

function [M, n] = normc(M)
   n = sqrt(sum(M.^2,1));     % Compute norms of columns
   M = bsxfun(@rdivide,M,n);  % Normalize M
   n = reshape(n,[],1);       % Store column vector of norms
end