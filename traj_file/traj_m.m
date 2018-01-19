%maneuver

if t==0
    %-------------- init --------------%
    T = 280;
    n = T/dt*2+1;
	angle = zeros(n,3); %deg
    speed = zeros(n,3); %m/s
    %---------------------------------------------------------------------%
    v0 = [1, 0, 90]; %m/s,deg, [horizontal velocity, down velocity, velocity direction]
    att0 = [90, 0, 0]; %deg, [psi,theta,gamma]
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
    if 60<t && t<=70
        angle(k,1) = 90*(70-t)/10;
    elseif 90<t && t<=100
        angle(k,1) = -90*(t-90)/10;
    elseif 120<t && t<=130
        angle(k,1) = -90*(130-t)/10;
    elseif 150<t && t<=160
        angle(k,1) = 90*(t-150)/10;
    elseif 180<t && t<=190
        angle(k,1) = 90*(190-t)/10;
    elseif 210<t && t<=220
        angle(k,1) = -90*(t-210)/10;
    end

    %-------------- pith --------------%
    angle(k,2) = 10*sin(pi*t/20)+40;
    
    %-------------- roll --------------%
    angle(k,3) = 50*sin(pi*t/60)+20;
    
    %-------------- vh --------------%
    
    %-------------- vd --------------%
    
    %-------------- vy --------------%
    if 60<t && t<=70
        vy = 90*(70-t)/10;
    elseif 90<t && t<=100
        vy = -90*(t-90)/10;
    elseif 120<t && t<=130
        vy = -90*(130-t)/10;
    elseif 150<t && t<=160
        vy = 90*(t-150)/10;
    elseif 180<t && t<=190
        vy = 90*(190-t)/10;
    elseif 210<t && t<=220
        vy = -90*(t-210)/10;
    end
    
    %*********************************************************************%
    speed(k,:) = [vh*cosd(vy), vh*sind(vy), vd];
end