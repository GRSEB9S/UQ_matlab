function beam
clc;
clear;
close all;

nRandomVarsUsed = 5; % number of random variables in K
t_final=2;           % time marching coeffs
dt=0.01;
T = 0:dt:t_final;
Method_TM2='Av';     % acceleration method for Newmark's scheme, Av(Average/FG(Fox&Goodwin)
N=100;               % Monte Carlo sample size
GS_ord = 3;
Kmodes = [10,15,20]; % number of modes used in GSD scheme
N_snap = 200;

load('geobeam.mat','M','K','C','F','nint','output','remain','remove','dfix');    %load M,K,C,F matrices

% Dynamic plot of the beam, assuming deterministic
[~,~] = FEM_det(T,Method_TM2);

% Plot the time-dependent part of the forcing
% Define the time variance of the applied force
time_dep = zeros(size(T));
a = 3/t_final;
b = 20*pi/t_final;
for i = 1:length(T)
    time_dep(i) = 1-exp(-a*(i-1)*dt)*(1-sin(b*(i-1)*dt));
end
plot(T,time_dep,'k-');
ylim([-1,1]);
xlabel('t(s)');
    
%initial conditions
u0 = zeros(nint,1);
ut0 = zeros(nint,1);

result = main(M,K,C,F,u0,ut0,t_final,dt,N_snap,nRandomVarsUsed,N,GS_ord,Kmodes,nint,remain,remove,dfix,Method_TM2,output);

save('beam.mat');

%% plot solutions
fprintf('-------------------------------------------------------------\n');
fprintf('Plotting...\n');
fprintf('-------------------------------------------------------------\n');
% load MCS results
load('MCS_beam_1e6.mat','mean_MC','std_MC');
%change mean_MC to result.mean_MCS if want to plot error against the
%freshly generated MCS results instead of the loaded MCS results

% Figure2/3: abs. err. of mean/std of gPC1/gPC2/gPC3
Yname{1}='gPC,p=1';
Yname{2}='gPC,p=2';
Yname{3}='gPC,p=3';
Linestyle{1}=':';
Linestyle{2}='-';
Linestyle{3}='--';
ploterr(3,T,result.mean_GS,mean_MC,Yname,'Absolute error in mean',Linestyle);
ploterr(3,T,result.std_GS,std_MC,Yname,'Absolute error in standard deviation',Linestyle);

% Figure4/5: abs. err. of mean/std. of gPC2/GSD10/GSD15/GSD20
Yname{1}='gPC,p=2';
Yname{2}='GSD,Kmodes=10';
Yname{3}='GSD,Kmodes=15';
Yname{4}='GSD,Kmodes=20';
Linestyle{1}='-';
Linestyle{2}=':';
Linestyle{3}='--';
Linestyle{4}='-.';
Ydata(1,:)=result.mean_GS(2,:);
Ydata(2:4,:)=result.mean_GSD(1:3,:);
ploterr(4,T,Ydata,mean_MC,Yname,'Absolute error in mean',Linestyle);
Ydata(1,:)=result.std_GS(2,:);
Ydata(2:4,:)=result.std_GSD(1:3,:);
ploterr(4,T,Ydata,std_MC,Yname,'Absolute error in standard deviation',Linestyle);

% Figure6/7: abs. err. of mean/std. of gPC2/AAPG1/AAPG2
Yname{1}='1st order AAPG';
Yname{2}='2nd order AAPG';
Yname{3}='gPC,p=2';
Linestyle{1}=':';
Linestyle{2}='--';
Linestyle{3}='-';
Ydata(1:2,:)=result.mean_AAPG(1:2,:);
Ydata(3,:)=result.mean_GS(2,:);
ploterr(3,T,Ydata,mean_MC,Yname,'Absolute error in mean',Linestyle);
Ydata(1:2,:)=result.std_AAPG(1:2,:);
Ydata(3,:)=result.std_GS(2,:);
ploterr(3,T,Ydata,std_MC,Yname,'Absolute error in standard deviation',Linestyle);

