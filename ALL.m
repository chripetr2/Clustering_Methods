
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

n_PC=3;%number of PC
 [eigenval,eigenvec,explain,X,mean_vec]=pca_fun(X',n_PC);
 %Initilization of the cluster representatives to start from a data point
max_m=8;%number of clusters

 [most_dist, theta_init] = most_dist_repre(X,max_m);
  X=X';

[px,nx]=size(X); % px= no. of rows (pixels) and nx=no. of columns (bands)

%%%  K-means     %%%%%
[theta,bel,J]=k_means(X',theta_init);

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

X_L=X';
subplot(4,3,1)
 plot3(X_L(1,bel==1),X_L(2,bel==1),X_L(3,bel==1),'.',...
X_L(1,bel==2),X_L(2,bel==2),X_L(3,bel==2),'*',X_L(1,bel==3),X_L(2,bel==3),X_L(3,bel==3),'o',...
X_L(1,bel==4),X_L(2,bel==4),X_L(3,bel==4),'x',X_L(1,bel==5),X_L(2,bel==5),X_L(3,bel==5),'d',...
X_L(1,bel==6),X_L(2,bel==6),X_L(3,bel==6),'p',X_L(1,bel==7),X_L(2,bel==7),X_L(3,bel==7),'ks',X_L(1,bel==8),X_L(2,bel==8),X_L(3,bel==8),'+')


subplot(4,3,2)
 imagesc(im_cl_label)

  Acc=[ tot_acc  ];
  c ={'K-means'};

subplot(4,3,3)
bar(  Acc)
text(1:length(Acc),Acc,num2str(Acc'),'vert','bottom','horiz','center'); 
box off
set(gca,'xticklabel',c)
ylabel('Accuracy')
ylim ([ 0 1   ])


%%%%% Fuzzy c-means %%%
q=2; %the fuzzyfier
[theta,U,obj_fun]=fuzzy_c_means(X',max_m,q);
[qw,bel2]=max(U');

n=150;
cl_label=bel2;
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


subplot(4,3,4)
 plot3(X_L(1,bel2==1),X_L(2,bel2==1),X_L(3,bel2==1),'.',...
X_L(1,bel2==2),X_L(2,bel2==2),X_L(3,bel2==2),'*',X_L(1,bel2==3),X_L(2,bel2==3),X_L(3,bel2==3),'o',...
X_L(1,bel2==4),X_L(2,bel2==4),X_L(3,bel2==4),'x',X_L(1,bel2==5),X_L(2,bel2==5),X_L(3,bel2==5),'d',...
X_L(1,bel2==6),X_L(2,bel2==6),X_L(3,bel2==6),'p',X_L(1,bel2==7),X_L(2,bel2==7),X_L(3,bel2==7),'ks',X_L(1,bel2==8),X_L(2,bel2==8),X_L(3,bel2==8),'+')


subplot(4,3,5)
 imagesc(im_cl_label)

  Acc=[ tot_acc  ];
  c ={'Fuzzy c-means'};

subplot(4,3,6)
bar(  Acc)
text(1:length(Acc),Acc,num2str(Acc'),'vert','bottom','horiz','center'); 
box off
set(gca,'xticklabel',c)
ylabel('Accuracy')
ylim ([ 0 1   ])

%%%%% PCMA %%%



 [theta,bel3]=pcma(X,theta_init',max_m);


n=150;
cl_label=bel3;
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


subplot(4,3,7)
 plot3(X_L(1,bel3==1),X_L(2,bel3==1),X_L(3,bel3==1),'.',...
X_L(1,bel3==2),X_L(2,bel3==2),X_L(3,bel3==2),'*',X_L(1,bel3==3),X_L(2,bel3==3),X_L(3,bel3==3),'o',...
X_L(1,bel3==4),X_L(2,bel3==4),X_L(3,bel3==4),'x',X_L(1,bel3==5),X_L(2,bel3==5),X_L(3,bel3==5),'d',...
X_L(1,bel3==6),X_L(2,bel3==6),X_L(3,bel3==6),'p',X_L(1,bel3==7),X_L(2,bel3==7),X_L(3,bel3==7),'ks',X_L(1,bel3==8),X_L(2,bel3==8),X_L(3,bel3==8),'+')


subplot(4,3,8)
 imagesc(im_cl_label)

  Acc=[ tot_acc  ];
  c ={'PCMA'};

subplot(4,3,9)
bar(  Acc)
text(1:length(Acc),Acc,num2str(Acc'),'vert','bottom','horiz','center'); 
box off
set(gca,'xticklabel',c)
ylabel('Accuracy')
ylim ([ 0 1   ])



%%% Possibilistic %%%%

 [ap,cp,mv,mc,iter,diffvec]=prob(X',theta_init,200);
 
theta=mv;
[~,bel4]=max(cp,[],2);

n=150;
cl_label=bel4;
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


subplot(4,3,10)
 plot3(X_L(1,bel4==1),X_L(2,bel4==1),X_L(3,bel4==1),'.',...
X_L(1,bel4==2),X_L(2,bel4==2),X_L(3,bel4==2),'*',X_L(1,bel4==3),X_L(2,bel4==3),X_L(3,bel4==3),'o',...
X_L(1,bel4==4),X_L(2,bel4==4),X_L(3,bel4==4),'x',X_L(1,bel4==5),X_L(2,bel4==5),X_L(3,bel4==5),'d',...
X_L(1,bel4==6),X_L(2,bel4==6),X_L(3,bel4==6),'p',X_L(1,bel4==7),X_L(2,bel4==7),X_L(3,bel4==7),'ks',X_L(1,bel4==8),X_L(2,bel4==8),X_L(3,bel4==8),'+')


subplot(4,3,11)
 imagesc(im_cl_label)

  Acc=[ tot_acc  ];
  c ={'Possibilistic'};

subplot(4,3,12)
bar(  Acc)
text(1:length(Acc),Acc,num2str(Acc'),'vert','bottom','horiz','center'); 
box off
set(gca,'xticklabel',c)
ylabel('Accuracy')
ylim ([ 0 1   ])

