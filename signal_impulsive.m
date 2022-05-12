function [X,S]=signal_impulsive(M, alpha, SNR, K,c)


N_alpha=length(alpha);
A=exp(-i*pi*(0:M-1)'*sin(alpha*pi/180));
%Vj=diag(sqrt((   (10*ones(N_alpha,1)).^(SNR/10)   )*  (100-99*(1-c))  /2   ));
Vj=diag(sqrt((   (10*ones(N_alpha,1)).^(SNR/10)   ) /2   ));
S=Vj*(randn(N_alpha,K)+i*randn(N_alpha,K));


Mu = [0;0];
Sigma = cat(3,[1],[100]);
P=[1-c,c];
gm = gmdistribution(Mu,Sigma,P);
r1=random(gm,M*K);
r2=random(gm,M*K);
r1=reshape(r1,M,K);
r2=reshape(r2,M,K);

noise=sqrt(1/2)*(r1+i*r2);
Z=noise;


X=A*S+noise;
