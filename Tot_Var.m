function Dx = Tot_Var(e,W)  
Dx = W;
global nonzero_1
global nonzero_2 
global nonzero_vector 
 
de = e(nonzero_1) - e(nonzero_2); 
Dx(nonzero_vector) = W(nonzero_vector).*de;

% for k = 1 : length(I)
%     i = I(k);
%     j = J(k);
%     Dx(j,i) = (e(j) - e(i))*W(i,j); 
% end
% clear de;
