function Out_upsample = Stochastic_upample(In, Indx_sample, M, N) 
    Out_upsample = zeros(M,N);
    Out_upsample(Indx_sample) =  In; 
end