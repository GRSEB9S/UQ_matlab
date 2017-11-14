function hex
clc;
clear;
close all;

nRandomVarsUsed = 5; % number of random variables in K
t_final=1;           % time marching coeffs
dt=0.005;
Method_TM2='Av';     % acceleration method for Newmark's scheme, Av(Average/FG(Fox&Goodwin)
N=100;               % Monte Carlo sample size
GS_ord = 2;
Kmodes = [10];       % number of modes used in GSD method
N_snap = 200;        % steps between snap shots

load('geohex5.mat','M','K','C','F','nint','output','remain','remove','dfix');

%initial conditions
u0 = zeros(nint,1);
ut0 = zeros(nint,1);

result = main(M,K,C,F,u0,ut0,t_final,dt,N_snap,nRandomVarsUsed,N,GS_ord,Kmodes,nint,remain,remove,dfix,Method_TM2,output);

save('hex.mat');

%% plot solutions
fprintf('-------------------------------------------------------------\n');
fprintf('Plotting...\n');
fprintf('-------------------------------------------------------------\n');
% load MCS results
load('MCS_hex_1e5.mat');%% plot solutions
%change mean_MC to result.mean_MCS if want to plot error against the
%freshly generated MCS results instead of the loaded MCS results
T = 0:dt:t_final;

% Figure9/10: abs. err. of mean/std. of AAPG1/AAPG2/gPC2/GSD10
Yname{1}='1st order AAPG';
Yname{2}='2nd order AAPG';
Yname{3}='gPC,p=2';
Yname{4}='GSD,Kmodes=10';
Linestyle{1}='-.';
Linestyle{2}='--';
Linestyle{3}='-';
Linestyle{4}=':';
Ydata(1:2,:)=result.mean_AAPG(1:2,:);
Ydata(3,:)=result.mean_GS(2,:);
Ydata(4,:)=result.mean_GSD(1,:);
ploterr(3,T,Ydata,mean_MC,Yname,'Absolute error in mean',Linestyle);
Ydata(1:2,:)=result.std_AAPG(1:2,:);
Ydata(3,:)=result.std_GS(2,:);
Ydata(4,:)=result.std_GSD(1,:);
ploterr(3,T,Ydata,std_MC,Yname,'Absolute error in standard deviation',Linestyle);




