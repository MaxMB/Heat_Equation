clear; close; clc;
%%% System constants
L = 1; % bar length
T = 10; % max time
a = 1; % thermal conductivity
K = [0.5,0.5]; % source constants
na = 0.5; % noise amplitude

%%% Simulation constants
m = 100; % # gaps x
h = L/m; % step x
x = 0:h:L; % spacial vector
n = 250; % # gaps t
k = T/n; % step t
t = 0:k:T; % time vector
ks = 0.05; % k step for Hypercube
Nmcm = 400; % # Monte Carlo estimations

%%% Finite difference method (FDM) solution
u = fdms(a,m,h,x,n,k,t,K);
un = u + na * rand(m+1,n+1) - na * rand(m+1,n+1); % uniformly distributed noise input signal [-na,+na]

%%% K vector estimation by FDM
Ke_hc = ls_hcm_fdm(un,a,m,h,x,n,k,t,ks); % hypercube
ue_hc = fdms(a,m,h,x,n,k,t,Ke_hc);
Ke_mc = ls_mcm_fdm(un,a,m,h,x,n,k,t,Nmcm); % Monte Carlo
ue_mc = fdms(a,m,h,x,n,k,t,Ke_mc);

%%% Plot output
subplot(2,2,1); mesh(t,x,u); title(['FDM analitical: k1=' num2str(K(1)) ' k2=' num2str(K(2))]);
    xlabel('t(s)'); ylabel('x(m)'); zlabel('T(ºC)');
subplot(2,2,2); mesh(t,x,un); title(['FDM noise amplitude = ' num2str(na)]);
    xlabel('t(s)'); ylabel('x(m)'); zlabel('T(ºC)');
subplot(2,2,3); mesh(t,x,ue_hc); title(['FDM hypercube: k1=' num2str(Ke_hc(1)) ' k2=' num2str(Ke_hc(2))]);
    xlabel('t(s)'); ylabel('x(m)'); zlabel('T(ºC)');
subplot(2,2,4); mesh(t,x,ue_mc); title(['FDM Monte Carlo: k1=' num2str(Ke_mc(1)) ' k2=' num2str(Ke_mc(2))]);
    xlabel('t(s)'); ylabel('x(m)'); zlabel('T(ºC)');
