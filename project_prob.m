

clear all
clc
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


%PCA functin the the data points
% X=normalize(X);

n_PC=2;%number of PC
 [eigenval,eigenvec,explain,X,mean_vec]=pca_fun(X',n_PC);
 %Initilization of the cluster representatives to start from a data point
max_m=8;
 [most_dist, theta_init] = most_dist_repre(X,max_m);
  X=X';

[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)


    
 [ap,cp,mv,mc,iter,diffvec]=prob(X',theta_init,50);
 
theta=mv;
[~,bel]=max(cp,[],2);

% The following code can be used after the execution of an algorithm 
% Let "label" be the px-dimensional vector, whose i-th element is the label
% of the class where the i-th vector in X has been assigned
% The code below, helps in depicting the results as an image again
% cl_label=(sum((rand(px,8)<ones(px,1)*[0 .1 .2 .3 .4 .5 .6 .7])'))'; %DO NOT pay attention in this line. It produces a fake labeling

cl_label=bel;
cl_label_tot=zeros(p*n,1);
cl_label_tot(existed_L)=cl_label;
im_cl_label=reshape(cl_label_tot,p,n);

%Compute accuracy

X_label=L(existed_L,:);


   x = unique(X_label);
   N = numel(x);
   count = zeros(N,1);
   for k = 1:N
      count(k) = sum(X_label==x(k));
   end
[Y,I]=sort(count,'descend');

for k=1:8
    m=I(k);
    tmp=[];

    for i=1:150
    
    for j=1:150
        
        if Salinas_Labels(i,j)==m
            
              tmp=[tmp;im_cl_label(i,j)];
              
        end    
    end
    end
    
    
     if k==1
       
        cl_name(k)=mode(tmp);
        acc(k)= sum(tmp==cl_name(k))/length(tmp);
        
   else
       
       for n=1:length(tmp)
           if sum(tmp(n)==cl_name(1:(k-1)))>0
            tmp(n)=0;
           end
       end
       s=nonzeros(tmp);
       cl_name(k)=mode(s);
       acc(k)= sum(tmp==cl_name(k))/length(tmp);
       
   end
    
end


tot_acc=sum(acc)/8;

subplot(1,2,1)
 imagesc(im_cl_label)

  Acc=[ tot_acc  ];
  c ={'PCMA'};

subplot(1,2,2)
bar(  Acc)
text(1:length(Acc),Acc,num2str(Acc'),'vert','bottom','horiz','center'); 
box off
set(gca,'xticklabel',c)
ylabel('Accuracy')
ylim ([ 0 1   ])

X=X';
figure(2), plot(X(1,bel==1),X(2,bel==1),'r.',...
X(1,bel==2),X(2,bel==2),'g*',X(1,bel==3),X(2,bel==3),'bo',...
X(1,bel==4),X(2,bel==4),'cx',X(1,bel==5),X(2,bel==5),'md',...
X(1,bel==6),X(2,bel==6),'yp',X(1,bel==7),X(2,bel==7),'ks',X(1,bel==8),X(2,bel==8),'*')
hold on
figure(2), plot(theta(1,:),theta(2,:),'k+')
figure(2)
hold on
% 
% X=X';
% figure(2), plot3(X(1,:),X(2,:),'k.')
% figure(2), plot3(X(1,bel==1),X(2,bel==1),X(3,bel==1),'.',...
% X(1,bel==2),X(2,bel==2),X(3,bel==2),'*',X(1,bel==3),X(2,bel==3),X(3,bel==3),'o',...
% X(1,bel==4),X(2,bel==4),X(3,bel==4),'x',X(1,bel==5),X(2,bel==5),X(3,bel==5),'d',...
% X(1,bel==6),X(2,bel==6),X(3,bel==6),'p',X(1,bel==7),X(2,bel==7),X(3,bel==7),'ks',X(1,bel==8),X(2,bel==8),X(3,bel==8),'+')
% hold on
