%bias+noise+walk

dts = dt*2;
n = floor((size(imu,1)-1)/2); %the number of inertial solving

sigma_psi = 0.1/3; %deg
sigma_v = 0.01/3; %m/s

%-------------------------------------------------------------------------%
Nx = 9; %state variables
Np = 9; %process noise
Nm = 4; %measure noise

X = zeros(Nx,1);
P = diag([[1,1,1]*1, [1,1,1]*1, [1,1,1]*1e-4]);
Q = diag([[1,1,1]*gyro_noise/dt, [1,1,1]*acc_noise/dt, [1,1,1]*3e-7]);
R = diag([(sigma_psi/180*pi)^2, [1,1,1]*sigma_v^2]);

Lp = Nx+Np;
Lm = Nx;
[gamma_p, Wmp, Wcp] = UKF_para(Lp, 0.01, 2, 0);
[gamma_m, Wmm, Wcm] = UKF_para(Lm, 0.01, 2, 0);

ChiZ = zeros(Nm,2*Lm+1);
%-------------------------------------------------------------------------%

nav = zeros(n,9);
filter = zeros(n,3);
output_P = zeros(n,Nx);
data = [(traj(:,7)+randn(2*n+1,1)*sigma_psi)/180*pi, traj(:,4:6)+randn(2*n+1,3)*sigma_v];

for k=1:n
    kj = 2*k+1;
    
	gyro0 = imu(kj-2, 1:3)' /180*pi;
    gyro1 = imu(kj-1, 1:3)' /180*pi;
    gyro2 = imu(kj  , 1:3)' /180*pi;
    acc0  = imu(kj-2, 4:6)';
    acc1  = imu(kj-1, 4:6)';
    acc2  = imu(kj  , 4:6)';
    
    %---------------------------------------------------------------------%
    Z = data(kj,:)';
    %--Predict--%
    X = [X;zeros(Np,1)]; %augmentation
    P = [P,zeros(Nx,Np);zeros(Np,Nx),Q];
    Psr = chol(P)';
    ChiX = [X, X*ones(1,Lp)+gamma_p*Psr, X*ones(1,Lp)-gamma_p*Psr]; %sigma point
    for p=1:(2*Lp+1)
        ChiX(:,p) = fun_state(ChiX(:,p), [gyro0;gyro1;gyro2;acc0;acc1;acc2], dts); %state equation
    end
    ChiX = ChiX(1:Nx,:);
    X = ChiX*Wmp;
    P = (ChiX-X*ones(1,2*Lp+1))*Wcp*(ChiX-X*ones(1,2*Lp+1))';
    %--Update--%
    Psr = chol(P)';
    ChiX = [X, X*ones(1,Lm)+gamma_m*Psr, X*ones(1,Lm)-gamma_m*Psr]; %sigma point
    for p=1:(2*Lm+1)
        ChiZ(:,p) = fun_measure(ChiX(:,p)); %measure equation
    end
    Zm = ChiZ*Wmm;
    Pxz = (ChiX-X*ones(1,2*Lm+1))*Wcm*(ChiZ-Zm*ones(1,2*Lm+1))';
    Pzz = (ChiZ-Zm*ones(1,2*Lm+1))*Wcm*(ChiZ-Zm*ones(1,2*Lm+1))' + R;
    K = Pxz/Pzz;
    X = X + K*(Z-Zm);
    P = P - K*Pzz*K';
    %---------------------------------------------------------------------%
    output_P(k,:) = sqrt(diag(P))';
    
    nav(k,4:6) = X(4:6)';
    nav(k,7:9) = X(1:3)' /pi*180;
    filter(k,1:3) = X(7:9)' /pi*180;
end

nav = [zeros(1,9); nav];

%*************************************************************************%
function [gamma, Wm, Wc] = UKF_para(L, alpha, beta, kappa)
    lamda = alpha^2*(L+kappa) - L;
    gamma = sqrt(L+lamda);
    Wm = ones(2*L+1,1)*0.5/(L+lamda);
    Wm(1) = lamda/(L+lamda);
    Wc = eye(2*L+1)*0.5/(L+lamda);
    Wc(1) = lamda/(L+lamda)+1-alpha^2+beta;
end

function x = fun_state(x, u, dt)
    gyro0 = u(1:3);
    gyro1 = u(4:6);
    gyro2 = u(7:9);
    acc0  = u(10:12);
    acc1  = u(13:15);
    acc2  = u(16:18);
    K1 = fun_dx(x, [gyro0;acc0]);
    K2 = fun_dx(x+dt*K1/2, [gyro1;acc1]);
    K3 = fun_dx(x+dt*K2/2, [gyro1;acc1]);
    K4 = fun_dx(x+dt*K3, [gyro2;acc2]);
    x = x + (K1+2*K2+2*K3+K4)*dt/6;
end

function z = fun_measure(x)
    z = [x(1);x(4:6)];
end

%x = [psi,theta,gamma, vx,vy,vz, ex,ey,ez, ngx,ngy,ngz, nax,nay,naz, nex,ney,nez]
%u = [wx,wy,wz, ax,ay,az]
function dx = fun_dx(x, u)
    da = [0, sin(x(3))/cos(x(2)), cos(x(3))/cos(x(2));
          0, cos(x(3)),          -sin(x(3));
          1, sin(x(3))*tan(x(2)), cos(x(3))*tan(x(2))] * (u(1:3)-x(7:9)-x(10:12));
    dv = angle2dcm(x(1),x(2),x(3))'*(u(4:6)-x(13:15)) + [0;0;10];
    dx = [da;dv;x(16:18);[0;0;0];[0;0;0];[0;0;0]];
end