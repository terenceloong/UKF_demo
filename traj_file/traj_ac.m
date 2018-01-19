%attitude change

if t==0
    %-------------- init --------------%
    T = 200;
    n = T/dt*2+1;
	angle = zeros(n,3); %deg
    speed = zeros(n,3); %m/s
    %---------------------------------------------------------------------%
    v0 = [0, 0, 0]; %m/s,deg, [horizontal velocity, down velocity, velocity direction]
    att0 = [0, 20, 20]; %deg, [psi,theta,gamma]
    %---------------------------------------------------------------------%
    Cnb = angle2dcm(att0(1)/180*pi, att0(2)/180*pi, att0(3)/180*pi);
    vh = v0(1);
    vd = v0(2);
    vy = v0(3);
    vn0 = [vh*cosd(vy), vh*sind(vy), vd];
    vb0 = vn0*Cnb';
    angle(1,:) = att0;
    speed(1,:) = vn0;
else
    angle(k,:) = angle(k-1,:);
    speed(k,:) = speed(k-1,:);
    
    %-------------- yaw --------------%

    %-------------- pith --------------%
    angle(k,2) = 10*sin(pi*t/20)+20;

    %-------------- roll --------------%
    angle(k,3) = 5*sin(pi*t/10)+20;

    %-------------- vh --------------%

    %-------------- vd --------------%

    %-------------- vy --------------%
    
    %*********************************************************************%
    speed(k,:) = [vh*cosd(vy), vh*sind(vy), vd];
end