function K = ls_mcm_fdm(un,a,m,h,x,n,k,t,Nmcm)
%%% Least squares by Monte Carlo method and finite difference method
% Initial estimation
K1=rand; K2=rand; % uniformly distributed draws
lsp = 0; % previous least square value
ue = fdms(a,m,h,x,n,k,t,[K1,K2]);
for i = 1:m+1
    for j = 1:n+1
        lsp = lsp + (un(i,j) - ue(i,j))^2;
    end
end
% Monte Carlo search estimation
for nn = 1:1:Nmcm-1
    k1=rand; k2=rand; % uniformly distributed draws in [0,1]
    ue = fdms(a,m,h,x,n,k,t,[k1,k2]);
    ls = 0; % least square
    for i = 2:m
        for j = 2:n+1
            ls = ls + (un(i,j) - ue(i,j))^2;
        end
    end
    if ls < lsp % found better estimation
        K1 = k1;
        K2 = k2;
        lsp = ls;
    end
end
K = [K1,K2];