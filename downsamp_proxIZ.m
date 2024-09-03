function e = downsamp_proxIZ(y,xti,xbark,M,N,Indx_sample,epsilon)
X = reshape(xti,M,N);
E = X; 
y_dowsamp = y - Stochastic_downsample(xbark, Indx_sample);  

up = Stochastic_upample(y_dowsamp, Indx_sample, M, N) ;  
E  = X+epsilon*up;

e = E(:); clear X x y E;
end
