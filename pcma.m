function [theta,bel]=pcma(X,theta_init,m)

[px,nx]=size(X);


eta=ones(1,m)*100+ randn(1,m)*1;%%0.7 and 0.01 for the normalized data 

iter_thresh=200;
iner_iter_thresh=100;
eU_thresh=0.05;
iter=0;
eU=100;
e_thresh=0.05;
%initialiazation of theta
theta=theta_init;

U=rand(px,m)*10;
rho=0.1;  
  

%PCMA
 while ((iter<iter_thresh) && (eU>eU_thresh))
    iter=iter+1
    U_old=U;
    dist_all=[];
    
    for j=1:m
        
          d(:,j)=sqrt( sum( ( ( (ones(px,1)*theta(j,:))'   -X' ).^2)));
        
         
          D(:,j)=abs(d(:,j));     
          
          U(:,j)=  exp(   -( D(:,j)/eta(j) ) );
   
    end
    
     
        in_iter=0;
        e=100;
    
        while (in_iter<iner_iter_thresh) && (e>e_thresh)
            
             d=d+(d==0)*10^(-10);
             theta_old=theta;
             in_iter = in_iter + 1;
             
             
             U1=(U .* (D ./ d))';
            
             
             
                  for j=1:m
                      
                          grad(j,:)= - U1(j,:)* (X - (ones(px,1)*theta(j,:)) );


                          theta(j,:) = theta(j,:)-rho*grad(j,:);
                          


                         d(:,j)=sqrt( sum( ( ( (ones(px,1)*theta(j,:))'   -X' ).^2)));

                         D(:,j)=abs(d(:,j));     

                         U(:,j)=  exp(   -( D(:,j)/eta(j) ) );
                        
                              %Computation on where the data points belong
        
%                            dist=sum(((ones(px,1)*theta(j,:))'   -X' ).^2);
%                          dist_all=[dist_all; dist];
                         
                  end
                                  
%                [q1,bel]=min(dist_all);  
         
         [~,bel]=max(U');
%          bel=bel';        
         e=sum(sum(abs(theta-theta_old)));
        end
    
    
    
    
    
   eU= sum(sum(abs(U-U_old)));
 end

end
    