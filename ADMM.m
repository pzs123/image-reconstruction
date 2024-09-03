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

% for i = 1:length(filePaths)
for i = 3:3
    I = imread(fullfile(folderTest,filePaths(i).name));
    XR = im2double(I(:,:,1));
    [M,N] = size(XR);
    xr = XR(:);
    xr = xr/max(xr); 
    X = reshape(xr,M,N);
    x = X(:);   

    x0 = ones(M*N,1);
    xbar0 = x0;
    z10 = ones(M*N,1);
    lambda1 = 0.02;
    lambda2 = 0.02;
    rep = 500;

    nu = 0.05;
    cu = 0.01;
    cz = 0.005;

    load(join(['.\Weight\W_weight_',num2str(M),'_',num2str(N),'.mat'])); %W

    A = Weig_A(W,M*N);
    segments = ['r-';'b-';'g-';'k-';'c-'];
    seg=segments(5,:);      

    Noise = rand(M*N,1);
    global SNR;
%     SNR = 15;
    var = sqrt(10^(-SNR/10)); 
    x_noise = x + var * Noise;

    XK = [];       
    errsal = []; 
    z20 = 0.02*ones(length(A(:,1)),1);
    pu = 0.02*ones(M*N,1);
    pz = 0.02*ones(length(A(:,1)),1);

    global rate
%     rate = 1.0;

    test_num = 1; 
    for K=1:test_num    

        RAndperm = randperm(M*N);
        Num_sample = round(M*N*rate);
        Indx_sample = RAndperm(1:Num_sample);
        y = Stochastic_downsample(x_noise, Indx_sample); %sampling,outorder
        clear x_noise
        
%         H = speye(M*N);
%         m = round(M*N*rate);
%         RAndperm = randperm(M*N);
%         Indx_sample = sort(RAndperm(1:m));
%         phi =  H(Indx_sample,:);
%         y = phi*x_noise;
%         invers = inv(phi'*phi+cu*H+cz*A'*A);  %Method 2: Inverse
        
        xk = x0;
        xbark = xbar0;
        z1k = z10;
        z2k = z20;

        errs = [];
        k = 1; 
        while k< rep

            xti = xk - nu*(-pu-cu*(z1k-xbark)-A'*pz-cz*A'*(z2k-A*xbark));
            xkb = xk;  
            xk = downsamp_proxIZ(y,xti,xbark,M,N, Indx_sample, nu);  %Gradient-Descent
            xbark = 2*xk-xkb;
            
%             xkb = xk; 
%             xk = invers*(phi'*y+pu+cu*z1k+A'*pz+cz*A'*z2k);
%             xbark = 2*xk-xkb;      %Method 2: Inverse


            As = A*xk;

            z1ti = xk-pu/cu;
            z2ti = As - pz/cz;
            for j = 1:M*N       %z1
                if abs(z1ti(j))-lambda1/cu>0
                    z1k(j,:) = sign(z1ti(j))*(abs(z1ti(j))-lambda1/cu);
                else
                    z1k(j,:) = 0;
                end
            end

            for j = 1:M*N       %z2
                if abs(z2ti(j))-lambda2/cz>0
                    z2k(j,:) = sign(z2ti(j))*(abs(z2ti(j))-lambda2/cz);
                else
                    z2k(j,:) = 0;
                end
            end

            pu = pu + cu*(z1k-xk);
            pz = pz + cz*(z2k-As);

            err = 10*log10(sum((xk-x).^2)/sum((x).^2));
            errs = [errs;err];

            k = k+1;
            k
        end
        XK = [XK,xk];   
        errsal = [errsal,errs]; 
        test_num
    end
    errsal = sum(errsal,2)/test_num;
    XK = sum(XK,2)/test_num;

    figure(10);
    Xf = reshape(XK,[M,N]);
    x = 1:k-1;
    hold on;
    plot(x,errsal,seg,'DisplayName','ADMM');
    ylabel('relative errors(db)');
    xlabel('Iterations');


    figure('NumberTitle', 'off', 'Name', 'ADMM');
    mat = Xf;
    imshow(mat,'border','tight','initialmagnification','fit');
    title(strcat('ADMM'));
    
%     imwrite(mat,join(['.\result\ADMM\samp\',num2str(i),'.jpg'])); 

end



