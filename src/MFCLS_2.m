function results = MFCLS_2(image,M)
    [x,y,z] = size(image);
    num_p = size(M,2);
    results = zeros(x,y,num_p);
    for i = 1:x
        for j = 1:y
            r = reshape(image(i,j,:),z,1);
            [ab] = MFCLS_pixel(r, M);
            ab = reshape(ab,[1,1,num_p]);
            results(i,j,:) = ab;
         end
    end
end

function a_fcls = MFCLS_pixel(r, M)
    p = size(M, 2);
    invMtM = inv(M'*M);
    a_ls = invMtM*M'*r;
    o = ones(p, 1);

    % step 1. set a_fcls = a_scls
    a_scls = a_ls + invMtM * o * inv(o' * invMtM * o) * (1 - o' * a_ls); %#ok<*MINV>
    a_fcls = a_scls;

    % get sign vectors
    s_ls = a_ls ./ abs(a_ls);
    s_scls = a_scls ./ abs(a_scls);

    % step 4. if any negative, then continue
    k = 0;
    while any(a_fcls < -0.00001) && k < 2
        % step 2. compute lambdas
        L = [1-o'*a_fcls; 1-s_ls'*a_fcls];
        R = [o'*invMtM*o o'*invMtM*s_ls; s_ls'*invMtM*o s_ls'*invMtM*s_ls];
        lambda = R \ L;
        
        % step 3. compute new a_fcls
        a_fcls = a_scls + invMtM * (lambda(1) * o + lambda(2) * s_scls);
        
        k = k + 1;
    end
end
