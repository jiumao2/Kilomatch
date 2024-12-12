function ind_out = findNearestPoint(seq, p)
% seq: 1xn double
% p: 1xm double

ind_out = zeros(size(p));
for ind_p = 1:length(p)
    p_this = p(ind_p);
    l = 1;
    r = length(seq);
    m = round((l+r)/2);
    
    flag = false;
    while r-l>1
        if seq(m) > p_this
            r = m;
        elseif seq(m) < p_this
            l = m;
        else
            ind_out(ind_p) = m;
            flag = true;
            break
        end
        m = round((l+r)/2);
    end

    if flag
        continue
    end
    
    if abs(seq(l)-p_this) < abs(seq(r)-p_this)
        ind_out(ind_p) = l;
    else
        ind_out(ind_p) = r;
    end
end
end