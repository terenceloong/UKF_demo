n = size(nav,1);
t = (0:n-1)*dts;

error = nav - traj(1:2:2*n,:);
for k=1:n
    if error(k,7)>300
        error(k,7) = error(k,7)-360;
    elseif error(k,7)<-300
        error(k,7) = error(k,7)+360;
    end
    if error(k,9)>300
        error(k,9) = error(k,9)-360;
    elseif error(k,9)<-300
        error(k,9) = error(k,9)+360;
    end
end

figure

subplot(3,2,1)
plot(t, error(:,4))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_x\rm(m/s)')
grid on

subplot(3,2,3)
plot(t, error(:,5))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_y\rm(m/s)')
grid on

subplot(3,2,5)
plot(t, error(:,6))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_z\rm(m/s)')
grid on

subplot(3,2,2)
plot(t, error(:,7))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\psi(\circ)')
grid on

subplot(3,2,4)
plot(t, error(:,8))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\theta(\circ)')
grid on

subplot(3,2,6)
plot(t, error(:,9))
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\gamma(\circ)')
grid on