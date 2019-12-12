clear
format compact
close all

load Salinas_Data

[p,n,l]=size(Salinas_Image); % Size of the Salinas cube


% Making a two dimensional array whose rows correspond to the pixels and
% the columns to the bands, containing only the pixels with nonzero label.

X_total=reshape(Salinas_Image, p*n,l);
X_total=abs(X_total);
L=reshape(Salinas_Labels,p*n,1);
existed_L=(L>0);   %This contains 1 in the positions corresponding to pixels with known class label
X=X_total(existed_L,:);
% X=normalize(X);

Y_total=reshape(Salinas_Labels, p*n,1);
Y_total=abs(Y_total);
Y=Y_total(existed_L,:);
bel=Y;
%PCA functin the the data points

m=3;

 [eigenval,eigenvec,explain,X,mean_vec]=pca_fun(X',m);
 X=X';
 
[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)





cl_label=bel;
cl_label_tot=zeros(p*n,1);
cl_label_tot(existed_L)=cl_label;
im_cl_label=reshape(cl_label_tot,p,n);


X=X';
% figure(2), plot3(X(1,:),X(2,:),'k.')
figure(2), plot3(X(1,bel==1),X(2,bel==1),X(3,bel==1),'.',...
X(1,bel==2),X(2,bel==2),X(3,bel==2),'*',X(1,bel==3),X(2,bel==3),X(3,bel==3),'o',...
X(1,bel==4),X(2,bel==4),X(3,bel==4),'x',X(1,bel==5),X(2,bel==5),X(3,bel==5),'d',...
X(1,bel==6),X(2,bel==6),X(3,bel==6),'p',X(1,bel==7),X(2,bel==7),X(3,bel==7),'ks',X(1,bel==8),X(2,bel==8),X(3,bel==8),'+')
hold on


% Creating the first three principal components
X=X-min(min(X)); % Making Y>=0
PCmat=zeros(3,p*n);
PCmat(1,existed_L)=X(1,:);
PCmat(2,existed_L)=X(2,:);
PCmat(3,existed_L)=X(3,:);
PCmat=PCmat';
PCmat=reshape(PCmat,p,n,3);
PCmat=(PCmat-min(min(min(PCmat))) )/( max(max(max(PCmat)))-min(min(min(PCmat))) );

% Depicting the first PC's of the salinas image cube
for i=1:3
    figure(i+11), imagesc(PCmat(:,:,i)); colorbar; axis off; axis image
    pause(1)
end