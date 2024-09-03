function z1k = z1_k(z1ti)
    z1k = z1ti;
    z1k(z1k>1) = 1;
    z1k(z1k<-1) = -1; 