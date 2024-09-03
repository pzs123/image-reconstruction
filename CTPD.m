clear   
% close all    

%%% read images
% folderTest  = '.\set14\gray';
folderTest  = '.\set5\color';
ext         =  {'*.jpg','*.png','*.bmp'};
filePaths   =  [];
for i = 1 : length(ext)
    filePaths = cat(1,filePaths, dir(fullfile(folderTest,ext{i})));
end

% for i = 1:length(filePaths)   %change
for i = 3:3                     %change
    I = imread(fullfile(folderTest,filePaths(i).name));
    XR = im2double(I(:,:,1));
    [M,N] = size(XR);
    xr = XR(:);
    xr = xr/max(xr); 
    X = reshape(xr,M,N);
    x = X(:);   %Convert to vector

    figure('NumberTitle', 'off', 'Name', 'Original image');
    mat = X;
    imshow(mat,'border','tight','initialmagnification','fit'); 

    global SNR;
%     SNR = 15;        %debug

    rand('seed',0);
    var = sqrt(10^(-SNR/10)); 
    Noise = rand(M*N,1);
%     Noise = randn(M*N,1);
    x_noise = x + var*Noise;  

    global x0;
    global xbar0;
    global z10;
    global Z20;
    global nu;
    global tau;
    global mu;
    global lambda1;
    global lambda2;
    global test_num;
    global rep;
    global Indx_sample;

    x0 = 1*ones(M*N,1);
    xbar0 = x0;
    z10 = 1*ones(M*N,1);
    Z20 = sparse(M*N,M*N);
    nu = 0.05;   
    mu = 0.1; 
    tau = mu;   %Step size
    sit = 1;   
    load(join(['.\Weight\W_weight_',num2str(M),'_',num2str(N),'.mat']));
%     W = Weight(X,sit);
%     save(join(['Weight\W_weight_',num2str(M),'_',num2str(N),'.mat']),'W');

    lambda1 = 0.02;
    lambda2 = 0.02;
    segments = ['r-';'b-';'g-';'k-';'c-'];
    seg=segments(1,:);      

    global nonzero_1
    global nonzero_2 
    global nonzero_vector 
    nonzero_vector = find(W~=0); 
    [nonzero_1,nonzero_2] = find(W~=0); 

    XK = [];       
    test_num = 1;   
    rep = 500;

    errsal = []; 

    global rate
%     rate = 1.0; %sampling,outorder

    for K=1:test_num 

        RAndperm = randperm(M*N);
        Num_sample = round(M*N*rate);
        Indx_sample = RAndperm(1:Num_sample);
        y = Stochastic_downsample(x_noise, Indx_sample); %sampling,outorder
        clear x_noise

        xk = x0;
        xbark = xbar0;
        z1k = z10;
        Z2k = Z20;

        errs = [];
        k = 1; 
        while k< rep

            divZ2 = div(W,Z2k);
            xti = xk - nu*(lambda1*z1k-lambda2*divZ2);
            xkb = xk;
            xk = downsamp_proxIZ(y,xti,xbark,M,N, Indx_sample, nu);  %Gradient-Descent
            xbark = 2*xk-xkb;

            dxk = Tot_Var(xbark,W);

            z1ti = z1k + tau * xbark;
            Z2ti = Z2k + mu * dxk;

            z1k = z1_k(z1ti);
            Z2k = Z2_k(Z2ti);

            err = 10*log10(sum((xk-x).^2)/sum((x).^2));
            errs = [errs;err];

            k = k+1;
            k
        end
        XK = [XK,xk];   
        errsal = [errsal,errs]; 
        K
    end
    errsal = sum(errsal,2)/test_num;
    XK = sum(XK,2)/test_num;

    figure(10);
    Xf = reshape(XK,[M,N]);
    x = 1:size(errsal,1);
    hold on;
    box on;
    plot(x,errsal,seg,'DisplayName','CTPD');
    axis normal;
    ylabel('relative errors(db)');
    xlabel('Iterations');
    legend('show');

    figure('NumberTitle', 'off', 'Name', 'CTPD');
    mat = Xf;
    imshow(mat,'border','tight','initialmagnification','fit'); 
    
%     imwrite(mat,join(['.\result\CTPD\samp\',num2str(i),'.jpg'])); 
end



