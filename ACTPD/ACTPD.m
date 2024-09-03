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

for J = 3:3
    I = imread(fullfile(folderTest,filePaths(J).name));
    XR = im2double(I(:,:,1));
    % XR = cell2mat(Train(21,1))'; %21
    [M,N] = size(XR);
    xr = XR(:);
    xr = xr/max(xr); 
    X = reshape(xr,M,N);
    x = X(:);   

    global Indx_sample;
    x0 = 1*ones(M*N,1);
    xbar0 = x0;
    z10 = 1*ones(M*N,1);
    Z20 = sparse(M*N,M*N);
    nu = 0.01;  %0.05
    mu = 0.1;   %0.1

  
    load(join(['.\Weight\W_weight_',num2str(M),'_',num2str(N),'.mat']));
    lambda1 = 0.02;  %0.02
    lambda2 = 0.02;  %0.02
    segments = ['-r';'-b';'-g';'-k';'-c'];
    seg=segments(1,:);      

    global nonzero_1
    global nonzero_2 
    global nonzero_vector 
    nonzero_vector = find(W~=0); 
    [nonzero_1,nonzero_2] = find(W~=0); 

    XK = [];       
    test_num = 1;   
    m=2;
%     rand('seed',0);
    Noise = rand(M*N,test_num); 

    rep = 500;
    rate = 1.00; 
    SNR = 15;
    for Tri=1:test_num  

        var = sqrt(10^(-SNR/10));  
        x_noise = x + var*Noise(:,Tri);

        RAndperm = randperm(M*N);
        Num_sample = round(M*N*rate);
        Indx_sample = RAndperm(1:Num_sample);
        y = Stochastic_downsample(x_noise, Indx_sample); %sampling,outorder
        clear x_noise

        xk(:,1) = x0; %求x1  
        z1k = z10;  %求z1k 
        Z2k = Z20;  %求Z2k 

    %     errs = 1; 
        errs = [];
        Rall = [];
        for k = 1: rep-1 
            mk = min(m,k);

            divZ2 = div(W,Z2k);

            xk_g(:,k) = g_x(xk(:,k),z1k,y(:,Tri),divZ2, nu, lambda1, lambda2,M,N);
            xbark = 2*xk_g(:,k)-xk(:,k);

            z1k_t = g_h(xbark,z1k, mu);
            if mk < k
                z1k_g = z1k_g(:,2:end);
                z1k_g(:,mk) = z1k_t;
            else
                z1k_g(:,mk) = z1k_t;
            end

            dxk = Tot_Var(xbark,W);
            Z2k_t = g_z(dxk,Z2k,mu);
            if mk < k && mk~=1
                Z2k_g = Z2k_g(2:end);
                Z2k_g{mk} = Z2k_t;
            else
                Z2k_g{mk} = Z2k_t;
            end

            R1 = xk_g(:,k)-xk(:,k);
            R2 = z1k_g(:,mk)-z1k;

            R3 = mean(Z2k_g{mk}-Z2k,2); 
            R0 = [R1,R2,R3];
            r = R0(:);

            Rall = [Rall,r];

            Rk  = Rall(:,k:-1:k-mk+1);
            Qkm = (Rk'*Rk)^(-1);
            one = ones(mk,1);
            ak = Qkm*one/(one'*Qkm*one);   

            if k<m
                for i = 1:k
                    aak1(:,i) = ak(i)*xk_g(:,k-i+1); 
                    aak2(:,i) = ak(i)*z1k_g(:,mk-i+1); 
                    aak3{i} = ak(i)*Z2k_g{mk-i+1}; 
                end
            else
                for i = 1:m
                    aak1(:,i) = ak(i)*xk_g(:,k-i+1); 
                    aak2(:,i) = ak(i)*z1k_g(:,mk-i+1); 
                    aak3{i} = ak(i)*Z2k_g{mk-i+1}; 
                end
            end


            xk(:,k+1) = sum(aak1,2);
            z1k = sum(aak2,2);
            Z2k = Z20;      %初始化求和
            for i = 1:length(aak3)
                Z2k = Z2k + aak3{i};
            end

            err = sum((xk(:,k+1)-x).^2)/sum(x.^2); 

            errs = [errs;err];
            k
        end
        XK = [XK,xk(:,end)];  
        errsal(Tri,:) = errs; 
        Tri
    end
    errsal = mean(errsal,1)/test_num; 
    XK = sum(XK,2)/test_num;

    figure(3); 
    x = 1:length(errsal);
    hold on;
    plot(x,errsal,seg);
    ylabel('RE(db)');
    xlabel('Iterations');

    Xf = reshape(XK,[M,N]);
    figure(4);
    mat = Xf;
    imshow(mat,'border','tight','initialmagnification','fit');  
%     
%     imwrite(mat,join(['.\result\ACTPD\A',num2str(J),'.jpg'])); 
end





