function [res]=Bayesian_DOA_Impulsive_joint(Y,search_area,N_alpha)

Y=Y/1000; % make a scale, where 1000 can be changed to other values.
[M,T]=size(Y);
K_hat=length(search_area);
reslu=search_area(2)-search_area(1);
search_mid_left=search_area-reslu/2;
search_mid_right=search_area+reslu/2;
Z=zeros(size(Y));
a_search=search_area*pi/180.;
A_theta=exp(-1i*pi*(0:M-1)'*sin(a_search));


%%%%%%%%%%%%%%
a=1 + 0.0001;
b=0.0001;
d=0.01;
maxiter=500;
tol=1e-3;
%initialization
beta=1;
delta=ones(K_hat,1);
gamma=ones(M,T);
mu=zeros(K_hat+M,T);
Sigma_all=zeros(K_hat+M,K_hat+M,T);
%%%%%%%%%%%%%%%
converged = false;
iter = 0;

while (~converged) || iter<=300 
       iter = iter + 1;
       delta_last = delta;
      %calculate mu and Sigma
       Phi=[A_theta, eye(M)];
       for ii=1:T
           deltai=[ delta; gamma(:,ii)];
           V_temp= 1/beta*eye(M) +  Phi*diag(deltai)*Phi';
           Vinv=inv(V_temp);
           Sigmai = diag(deltai) -diag(deltai) * Phi' * Vinv * Phi *  diag(deltai);
           mu(:,ii) = beta * Sigmai * Phi' * Y(:,ii);
           Sigma_all(:,:,ii)=Sigmai;
       end

       %update delta
       Temp1=0;
       for ii=1:T
           Temp1=Temp1+ Sigma_all(:,:,ii);   
       end
       temp=sum( mu.*conj(mu), 2) + real(diag(Temp1));     
       delta= ( -T+ sqrt(  T^2 + 4*d* real(temp(1:K_hat) ) ) ) / (  2*d   );     

       for ii=1:M
           for tt=1:T;
               Temp1= abs(mu(K_hat+ii,tt))^2  +  Sigma_all(K_hat+ii,K_hat+ii,tt); 
               gamma(ii,tt)=abs(-1 + sqrt( 1+4*d*Temp1  )  )  / (  2*d   );           
           end
       end

       %update beta
       temp=0;
       for ii=1:T
          temp=temp+sum(diag(Phi*Sigma_all(:,:,ii)*Phi'));
       end

       resid=Y-Phi*mu;
       beta=abs ( ( T*M + (a-1))/( b +  norm(resid, 'fro')^2  +  temp     ))  ;

        % stopping criteria
          erro=norm(delta - delta_last)/norm(delta_last);
        if erro < tol || iter >= maxiter
            converged = true;
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%off grid
        etc=M;
        mu2=mu(1:K_hat,:);
        f=sqrt(sum(mu2.*conj(mu2),2));
        [~,sort_ind]=sort(f);
        index_amp=sort_ind(end:-1:end-etc+1);
        for j=1:length(index_amp)            
                ii=index_amp(j);
                mut=mu(ii,:);
                Sigmat=Sigma_all(:,ii,:);
                phi=mut*mut' +  sum( Sigmat(ii,:,:));
                tempind=[1:K_hat+M]; tempind(ii)=[];
                Yti=Y- Phi(:, tempind)*mu(tempind,:);

                zer=0;
                for jj=1:T
                   zer=zer+  Phi(:,tempind)*  (Sigmat(tempind,:,jj));
                end
                varphi= zer  -Yti*(mut');
                z1=[1:M-1]';
                c=zeros(M,1);
                c(1)=M*(M-1)/2*phi;
                c(2:end)=z1.*varphi(2:end);            
                %%%%%%%%% root method 
               ro=roots(c);
               abs_root=abs(ro);
               [~,indmin]=min(abs(abs_root-1));
               angle_cand=  asin(-angle(ro(indmin))/pi)   /pi*180;      

               if angle_cand<=search_mid_right(ii) && angle_cand>=search_mid_left(ii)
                  search_area(ii)=angle_cand;
                  A_theta(:,ii)= exp(-1i*pi*(0:M-1)'*sin(angle_cand*pi/180));               
               end
         end        
end
mu=mu(1:K_hat,:);
Pm=mean(mu.*conj(mu),2);
Pm=Pm/max(Pm);
[res]=findmax_new(search_area,Pm,N_alpha,1);

