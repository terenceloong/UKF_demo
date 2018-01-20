%a simple UKF

T = 10;
dt = 0.01;

n = T/dt;
t = (1:n)'*dt;
sigma_z = 0.05;
data = sin(1.5*t) + randn(n,1)*sigma_z; %measure
filter = zeros(n,2);
output_P = zeros(n,2);

%-------------------------------------------------------------------------%
X = [0;0];
P = diag([1,1]);
R = sigma_z^2;
L = 2;
[gamma, Wm, Wc] = UKF_para(2, 0.01, 2, 0);
%-------------------------------------------------------------------------%

for k=1:n
    t = k*dt;
    
    %---------------------------------------------------------------------%
    Z = data(k);
    %--Predict--%
    Psr = gamma*chol(P)';
    ChiX = [X, X*ones(1,L)+Psr, X*ones(1,L)-Psr]; %sigma point
    for p=1:(2*L+1)
        ChiX(:,p) = fun_state(ChiX(:,p), t-dt, dt); %state equation
    end
    X = ChiX*Wm;
    P = (ChiX-X*ones(1,2*L+1))*Wc*(ChiX-X*ones(1,2*L+1))';
    %--Update--%
    Psr = gamma*chol(P)';
    ChiX = [X, X*ones(1,L)+Psr, X*ones(1,L)-Psr]; %sigma point
    ChiZ = [1,0]*ChiX; %measure equation
    Zm = ChiZ*Wm;
    Pxz = (ChiX- X*ones(1,2*L+1))*Wc*(ChiZ-Zm*ones(1,2*L+1))';
    Pzz = (ChiZ-Zm*ones(1,2*L+1))*Wc*(ChiZ-Zm*ones(1,2*L+1))' + R;
    K = Pxz/Pzz;
    X = X + K*(Z-Zm);
    P = P - K*Pzz*K';
    %---------------------------------------------------------------------%
    
    filter(k,:) = X';
    output_P(k,:) = sqrt(diag(P))';
end

t = (1:n)*dt;
figure
plot(t,data, 'LineWidth',1)
hold on
plot(t,filter(:,1), 'LineWidth',1)
grid on
legend('measure','after filtering')
xlabel('\itt\rm(t)')

figure
plot(t,filter(:,2), 'LineWidth',1)
grid on
xlabel('\itt\rm(t)')
title('\omega estimation')

%*************************************************************************%
function [gamma, Wm, Wc] = UKF_para(L, alpha, beta, kappa)
    lamda = alpha^2*(L+kappa) - L;
    gamma = sqrt(L+lamda);
    Wm = ones(2*L+1,1)*0.5/(L+lamda);
    Wm(1) = lamda/(L+lamda);
    Wc = eye(2*L+1)*0.5/(L+lamda);
    Wc(1) = lamda/(L+lamda)+1-alpha^2+beta;
end

function x = fun_state(x, t, dt)
    K1 = fun_dx(x, t);
    K2 = fun_dx(x+dt*K1/2, t+dt/2);
    K3 = fun_dx(x+dt*K2/2, t+dt/2);
    K4 = fun_dx(x+dt*K3, t+dt);
    x = x + (K1+2*K2+2*K3+K4)*dt/6;
end

function dx = fun_dx(x, t)
    dx = [x(2)*cos(x(2)*t);0];
end