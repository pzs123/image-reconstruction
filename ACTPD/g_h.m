function hkp = g_h(xbark,hk,mu)
    z1ti = hk + mu * xbark;
    hkp = z1_k(z1ti);