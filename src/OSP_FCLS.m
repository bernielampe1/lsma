function results = OSP_FCLS(HIM, M)
    p = size(M, 2);
    [xx, yy, ll] = size(HIM);

    results = zeros(xx, yy, p);

    projs = {};
    denoms = {};
    for i = 1:p
        m = M;
        m(:,i) = [];
        projs{end+1} = eye(ll) - m*inv(m'*m)*m';
    end

    for x = 1:xx
        for y = 1:yy
            r = squeeze(HIM(x, y, :));
            a = zeros(p, 1);
            for i = 1:p
                a(i) = norm(projs{i} * r);
            end
            results(x, y, :) = a./norm(a);
        end
    end 
end
