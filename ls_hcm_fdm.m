function K = ls_hcm_fdm(un,a,m,h,x,n,k,t,ks)
%%% Least squares by hypercube method and finite difference method
% Initial estimation
K1=rand; K2=rand;
lsp = 0; % previous least square value
ue = fdms(a,m,h,x,n,k,t,[K1,K2]);
for i = 1:m+1
    for j = 1:n+1
        lsp = lsp + (un(i,j) - ue(i,j))^2;
    end
end
% Hypercube search estimation
kf = 1; % final k value
ki = 0; % initial k value
for k1 = ki:ks:kf
    for k2 = ki:ks:kf
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
end
K = [K1,K2];