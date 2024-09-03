clear

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
    randn('state',4)
    
    ori = x;

    global rate;   
%     rate = 1.0; %sampling,outorder
    I = speye(M*N);
    RAndperm = randperm(M*N);
    Num_sample = round(M*N*rate);
    r1 = sort(RAndperm(1:Num_sample)); 
    A = I(r1,:);
    
    global SNR; 
%     SNR = 15;
    var = sqrt(10^(-SNR/10));   % var = 0.05; 
    rand('seed',0);
    Noise = rand(size(A,1),1);
    y = A*x + var*Noise;

    epsilon = 1e-10;
    itermax = 500;
    lambda = 0.02;
    x0 = ones(M*N,1);

    [x_2,error,errs] = cs_fista(y,A,lambda,epsilon,itermax,x);

    figure(10);
    hold on;
    box on;
    plot(errs,'m-','DisplayName','FISTA');
    axis normal;
    ylabel('relative errors(db)');
    xlabel('Iterations');
    
    xf = x_2;
    Xf = reshape(xf,M,N);
    figure('NumberTitle', 'off', 'Name', 'FISTA');
    mat = Xf;
    imshow(mat,'border','tight','initialmagnification','fit');  
    
%     imwrite(mat,join(['.\result\FISTA\samp\',num2str(i),'.jpg']));
    
end


