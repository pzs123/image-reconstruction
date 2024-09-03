function Z2k = Z2_k(Z2ti)
    Z2k = Z2ti;
    V_norm = vecnorm(Z2ti);
    logical_indices = V_norm > 1;
    Z2k(:,logical_indices) = Z2ti(:,logical_indices)./ V_norm(logical_indices);
    end
    
%     Z2k  = Z2ti;
%     s = sum(Z2ti.^2,2);
%     p = find(s>1);
%     for i = 1:length(p)
%         ii = p(i);
%         Z2k(ii,:) = Z2ti(ii,:)/sqrt(s(ii));
%     end
