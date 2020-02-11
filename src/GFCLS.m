function results = GFCLS(HIM, M)
    p = size(M, 2);
    [xx, yy, ~] = size(HIM);

    results = zeros(xx, yy, p);

    E = [ones(1, p);M];
    denom = abs(sqrt(det(E'*E)));

    for x = 1:xx
        for y = 1:yy
            r = squeeze(HIM(x, y, :));
            a = ratios(r, p, M) ./ denom;
            results(x, y, :) = a;
        end
    end 
end

function v = ratios(r, p, M)
    v = zeros(p, 1);
    for i = 1:p
        m = M;
        m(:,i) = r;
        v(i) = abs(sqrt(det(m'*m)));
    end
end
