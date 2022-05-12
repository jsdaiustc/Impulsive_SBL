clear;
close all;

alpha=[-30,0]+rand(1,2)*10;
N_alpha=length(alpha);
SNR=20;
K=20;
M=8;
resolution=2;
c=0.1;
[X,P_noise]=signal_impulsive(M, alpha, SNR, K,c);
search_area=[-90:resolution:90];
[res_joint]=Bayesian_DOA_Impulsive_joint(X,search_area,N_alpha);     % corresponds to Bayes-optimal method in paper. 
res_joint-alpha'
