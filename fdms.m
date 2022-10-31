function u = fdms(a,m,h,x,n,k,t,K) % Finite difference method solution
la = a^2 * k / h^2; % lambda
u = zeros(m+1,n+1); % FDM solution declaration
u(:,1) = ic(x); % initial condition

% A*w = C | A -> D (main diagonal) & R (ratio)
R = zeros(m,1); % ratio for LSS
D = zeros(m+1,1); % main diagonal
R(1) = 0.5 * la / (la+1);
D(1) = la + 1;
for i=2:m
    D(i) = la + 1 - 0.5 * la * R(i-1);
    R(i) = 0.5 * la / D(i);
end
D(m+1) = la + 1 - 0.5 * la * R(m);

% Linear system solver A*w(j+1) = B*w(j) + S = C
for j=2:n+1
    s = 0.5 * k * (source(t(j-1),K) + source(t(j),K));
    C = zeros(m+1,1);
    C(1) = (1+la)*u(1,j-1) - 0.5*la*u(2,j-1) + s;
    for i=2:m
        C(i) = (1+la)*u(i,j-1) - 0.5*la*(u(i-1,j-1)+u(i+1,j-1)) + s + R(i-1)*C(i-1);
    end
    C(m+1) = (1+la)*u(m+1,j-1) - 0.5*la*u(m,j-1) + s + R(m)*C(m);
    for i=m:-1:2 % Boundary conditions: w(t,0) = w(t,L) = 0
        u(i,j) = (C(i) + u(i+1,j)*0.5*la) / D(i);
    end
end