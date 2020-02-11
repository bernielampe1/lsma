function results = MFCLS(image,M,tol)
    [x,y,z] = size(image);
    num_p = size(M,2);
    results = zeros(x,y,num_p);
    for i = 1:x
        for j = 1:y
            r = reshape(image(i,j,:),z,1);

            % step 1. init
            delta = 1/(max(max(M)));        
            s = [delta.*r;1];
            N = [delta.*M;ones(1,num_p)];

            [ab] = MFCLS_pixel(s, N, tol);
            ab = reshape(ab,[1,1,num_p]);
            results(i,j,:) = ab;
         end
    end
end

function a_fcls = MFCLS_pixel(r, M, tol)
    %#ok<*MINV>
    % init
    p = size(M,2);
    k = 0;
    P = ones(1, p);
    R = zeros(1, p);

    % step 2.
    invMtM = inv(M'*M);
    a_ls = invMtM*M'*r;
    one = ones(p,1);
    a_scls = a_ls + invMtM * one * inv(one' * invMtM * one) * (1 - one' * a_ls);
    a_fcls = a_scls;

    % step 3. check for termination
    while abs(sum(abs(a_fcls)) - 1) > tol && k < 10
        % step 4.
        k = k + 1;

        % step 5. move indicies from P to R for negative a_fcls(i)
        for i = 1:p
            if a_fcls(i) < 0 && P(i) == 1
                P(i) = 0;
                R(i) = 1;
            end
        end

        c = 0;
        while c < 10
            c = c + 1;
            
            % step 6.
            a_R = a_scls(R == 1);
            
            % step 7.
            Phi_P = invMtM(R == 1, R == 1);
            Phi_R_c = invMtM(:, R == 1)' * ones(p, 1);
            Phi_R_r = invMtM * ones(p,1);
            Phi_R_r = Phi_R_r(R == 1);
            Phi_corner = ones(1,p)*invMtM*ones(p,1);
        
            Phi = [Phi_P Phi_R_c; Phi_R_r' Phi_corner]; % ambiguity about which corner

            % step 8
            lambda = Phi \ [a_R; 0];

            % step 9
            z = length(lambda);
            if all(lambda(1:z-1) < 0)
                break
            end

            % step 10
            [~, ind] = max(lambda);
            inds_R = find(R);
            ind = inds_R(ind);
            R(ind) = 0;
            P(ind) = 1;
        end

        % step 11
        Psi = invMtM(:, R == 1);

        % step 12
        l = length(lambda);
        lambda_1 = lambda(l);
        lambda_2 = lambda(1:l-1);

        % step 13
        a_fcls = a_scls - ones(1,p) * invMtM * ones(p, 1) * lambda_1 * Psi * lambda_2;
    end
    
    a_fcls(a_fcls < 0) = 0;
end
