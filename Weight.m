function W = Weight(X,sit) 
    [M,N] = size(X);
%     K = M*N;
%     for k1 = 1:K
%         for k2 = 1:K
%             i1 = mod(k1,M);
%             j1 = ceil(k1/M);
% 
%             i2 = mod(k2,M);
%             j2 = ceil(k2/M);
%             if i1==0
%                 i1=M;
%             end
%             if i2==0
%                 i2=M;
%             end
%             r = (i1-i2)^2+(j1-j2)^2;
%             W(k1,k2) = exp(-r/(2*sit^2));
%         end
%     end
%     W(W<0.001) = 0;
W = sparse(M*N,M*N);
w1 = exp(-1/(2*sit^2));     %相邻(距离平方)
w2 = exp(-2/(2*sit^2));     %对角
for m=1:M*N
    if m==1     %第一行第一列
        W(m,m+1) = w1; W(m,m+M) = w1;
        W(m,m+M+1) = w2;
    elseif m==M*N-M+1   %第一行最后一列
        W(m,m+1) = w1; W(m,m-M) = w1;
        W(m,m-M+1) = w2;
    elseif m==M     %最后一行第一列
        W(m,m-1) = w1; W(m,m+M) = w1;
        W(m,m+M-1) = w2;
    elseif m== M*N  %最后一行最后一列
        W(m,m-1) = w1; W(m,m-M) = w1;
        W(m,m-M-1) = w2;
    elseif mod(m,M)==1      %第一行
        W(m,m+1) = w1; W(m,m-M) = w1; W(m,m+M) = w1;
        W(m,m-M+1) = w2; W(m,m+M+1) = w2;
    elseif mod(m,M)==0   %最后一行
        W(m,m-1) = w1; W(m,m-M) = w1; W(m,m+M) = w1;
        W(m,m-M-1) = w2; W(m,m+M-1) = w2;
    elseif m<M   %第一列
        W(m,m-1) = w1; W(m,m+1) = w1; W(m,m+M) = w1;
        W(m,m+M-1) = w2; W(m,m+M+1) = w2;
    elseif m>M*N-M
        W(m,m-1) = w1; W(m,m+1) = w1; W(m,m-M) = w1;
        W(m,m-M-1) = w2; W(m,m-M+1) = w2; 
    else 
        W(m,m-1) = w1; W(m,m+1) = w1;
        W(m,m-M) = w1; W(m,m+M) = w1;
        W(m,m-M-1) = w2; W(m,m-M+1) = w2;
        W(m,m+M-1) = w2; W(m,m+M+1) = w2;
    end
end
