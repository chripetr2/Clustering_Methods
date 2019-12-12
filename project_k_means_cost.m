
clear
format compact
close all

load Salinas_Data

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube

%Depicting the bands of the Salinas cube
% for i=1:l
%     figure(1); imagesc(Salinas_Image(:,:,i))
%     pause(0.1)
% end

% Making a two dimensional array whose rows correspond to the pixels and
% the columns to the bands, containing only the pixels with nonzero label.

X_total=reshape(Salinas_Image, p*n,l);
X_total=abs(X_total);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to pixels with known class label
X=X_total(existed_L,:);
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)


%Initilization of the cluster representatives to start from a data point
% max_m=8;
for max_m=1:10

    for m=1:max_m

        k=randi([1 px]);

       theta_init(m,:)=X(k,:);

    end


    [theta,bel,J]=k_means(X',theta_init');

    [eigenval,eigenvec,explain,Y,mean_vec]=pca_fun(X',m);

% The following code can be used after the execution of an algorithm 
% Let "label" be the px-dimensional vector, whose i-th element is the label
% of the class where the i-th vector in X has been assigned
% The code below, helps in depicting the results as an image again
% cl_label=(sum((rand(px,8)<ones(px,1)*[0 .1 .2 .3 .4 .5 .6 .7])'))'; %DO NOT pay attention in this line. It produces a fake labeling

% cl_label=bel;
% cl_label_tot=zeros(p*n,1);
% cl_label_tot(existed_L)=cl_label;
% im_cl_label=reshape(cl_label_tot,p,n);
% figure(10), imagesc(im_cl_label)

    cost(max_m)=J;
end

mm=[1:10];
plot(mm,cost)