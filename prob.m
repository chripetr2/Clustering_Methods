function [ap,cp,mv,mc,iter,diffvec]=prob(X,theta_ini,maxiter)



    mv=theta_ini;

    [l,N]=size(X);
    [l,m]=size(mv);
    e = 0.01;
   
    
    mc=rand(l,l,m);
for i=1:m % initialize covariance matrix
   mc(:,:,i)= cov(X(:,randi([1 10]))');
end

    % Initialization of the conditional probabilities
    for i=1:N
        cp(i,:)=rand(1,m);
        sumcp=sum(cp(i,:));
        cp(i,:)=cp(i,:)/sumcp;
    end

    % GPrAS scheme - normal pdfs
    iter=0;
    diff=e+1;
    mvold=mv;            
    mcold=mc;            
    apold=ones(1,m)/m;   
    diffvec=[];
    while (iter<maxiter)&(diff>e)
        iter=iter+1
        diffvec=[diffvec diff];

        % a priori probabilities
        temp=sum(cp);
        ap=temp/N;

        % mean vectors
        for i=1:m
            tot=zeros(l,1);
            for j=1:N
                tot=tot+cp(j,i)*X(:,j);
            end
            mv(:,i)=tot/temp(i);
        end

        % covariance matrices
        for i=1:m
            tot=zeros(l);
            for j=1:N
                tot=tot+cp(j,i)*(X(:,j)-mv(:,i))*(X(:,j)-mv(:,i))';
            end
            mc(:,:,i)=tot/temp(i);
        end

        diff=sum(abs(ap-apold))+sum(sum(abs(mv-mvold)))+sum(sum(sum(abs(mc-mcold))));
        apold=ap;
        mvold=mv;
        mcold=mc;

        % the determinants of the covariance matrices
        for i=1:m
            dete(i)=det(mc(:,:,i));
        end

        % conditional probabilities
        for j=1:N
            for i=1:m
                gama(j,i)=dete(i)^(-.5)*exp(-.5*(X(:,j)-mv(:,i))'*inv(mc(:,:,i))*(X(:,j)-mv(:,i)) )*ap(i);
            end
            sumtot=sum(gama');
        end

        for j=1:N
           for i=1:m
               cp(j,i)=gama(j,i)/sumtot(j);
           end
        end
    end