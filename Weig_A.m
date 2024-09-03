function A = Weig_A(W,N)
    A = -W;
    for i=1:N

        A(i,i)=sum(W(i,:));

    end
    
% C=A;