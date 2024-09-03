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

    global SNR;
%     SNR = 15;
    rand('seed',0);
    var = sqrt(10^(-SNR/10));   % var = 0.05; 
    Noise = rand(M*N,1);
    x_noise = x + var*Noise;
    
    x0 = 1*ones(M*N,1);
    xbar0 = x0;
    z10 = 1*ones(M*N,1);
    global tau;
    global nu;
    global lambda1;
    global test_num;
    global rep;
    global Indx_sample;

    segments = ['r-';'b-';'g-';'k-';'c-'];
    seg=segments(3,:);      

    XK = [];       
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

        errs = [];
        k = 1; 
%         rep = 500;
        while k< rep
            xti = xk - nu*(lambda1*z1k);
            xkb = xk;   
            xk = downsamp_proxIZ(y,xti,xbark,M,N, Indx_sample, nu); 
            xbark = 2*xk-xkb;

            z1ti = z1k + tau * xbark; 

            z1k = z1_k(z1ti);

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
    plot(x,errsal,seg,'DisplayName','LASSO');
    ylabel('relative errors(db)');
    xlabel('Iterations');

    figure('NumberTitle', 'off', 'Name', 'LASSO');
    mat = Xf;
    imshow(mat,'border','tight','initialmagnification','fit');
    title(strcat('LASSO'));

%     imwrite(mat,join(['.\result\LASSO\samp\',num2str(i),'.jpg'])); 

end

