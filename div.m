function D = div(W,H) 
% [size_a,size_b] = size(W); 
% WH = sparse(size_a,size_b);
% global nonzero_vector 
% WH(nonzero_vector) = W(nonzero_vector).*H(nonzero_vector); 

WH     = W.*H; 
WH1    = sum(WH,1);
WH2    = sum(WH,2);
D      = WH1(:) - WH2(:);  
D      = sparse(D); 
clear WH WH1 WH2 ;
end
 