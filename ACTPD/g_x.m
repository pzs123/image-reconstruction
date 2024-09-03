function xkp = g_x(xk,hk,y,divZ2,nu, lambda1, lambda2,M,N)
    global Indx_sample;
    xti = xk - nu*(lambda1*hk-lambda2*divZ2); 
    xkp = downsamp_proxIZ(y,xti,xk,M,N, Indx_sample, nu);  %Gradient-Descent
 