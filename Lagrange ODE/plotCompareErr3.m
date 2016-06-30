function plotCompareErr3(filename)
% Plot results from compareErr3
% plot 1: objective error
% plot 2: trajectory error
% plot 3: max s error
% plot 4: max θ error
% x-axis: iteration
% legend: interior-point, sqp

if nargin >= 1
    load(filename);
    load('data/traj_x_comp.mat','x_comp','p','s','sd','sdd','theta','thetad','thetadd','t');
    load('data/traj_x_p0.mat','x_p0','p_0','s','sd','sdd','theta','thetad','thetadd','t');
    load('data/traj_x_ref.mat','x_ref','p_opt','s','sd','sdd','theta','thetad','thetadd','t');
else
    load 'data/compareErr3.mat';
end

% object names:
% absErr, relErr, absSTErr, relSTErr, fctVal, p_0, maxsdiff, maxtdiff

% font size
fsz = 15;

% Phi error
figure;
%subplot(2,2,1);
semilogy(squeeze(fctVal(1,1,:)),'-k');
hold on;
semilogy(squeeze(fctVal(2,1,:)),'-.k');
grid on;
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error', 'FontSize',fsz);
title('$0.5 \| \Phi (\bar{x},u,p) \|^2$','interpreter','latex', 'FontSize',fsz);

% x error
figure;
semilogy(squeeze(absErr(1,1,:)),'-k');
hold on;
semilogy(squeeze(absErr(2,1,:)),'-.k');
grid on;
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error', 'FontSize',fsz);
title('$0.5 \| \bar{x}-x(p)\|^2$','interpreter','latex', 'FontSize',fsz);

% max s error
figure;
semilogy(squeeze(maxsdiff(1,1,:)),'-k');
hold on;
semilogy(squeeze(maxsdiff(2,1,:)),'-.k');
grid on;
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error [m]', 'FontSize',fsz);
title('$\max_n(|\bar{s}_n - s_n(p)|)$','interpreter','latex', 'FontSize',fsz);

% max θ error
figure;
semilogy(squeeze(maxtdiff(1,1,:)),'-k');
hold on;
semilogy(squeeze(maxtdiff(2,1,:)),'-.k');
grid on;
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error [rad]', 'FontSize',fsz);
title('$\max_n(|\bar{\theta}_n - \theta_n(p)|)$','interpreter','latex', 'FontSize',fsz);

% s trajectories
figure;
plot(t,x_ref(1,:),'-k');
hold on;
plot(t,x_p0(1,:),'-.k');
plot(t,x_comp(1,:),'--k');
grid on;
h = legend('$\bar{s}$','$s(p_0)$','$s(p_{opt})$');
set(h,'Interpreter','latex','FontSize',fsz);
xlabel('time', 'FontSize',fsz);
ylabel('val [m]', 'FontSize',fsz);
title('$s$ for different parameters','interpreter','latex','FontSize',fsz);

% θ trajectories
figure;
plot(t,x_ref(2,:),'-k');
hold on;
plot(t,x_comp(2,:),'--k');
plot(t,x_p0(2,:),'-.k');
grid on;
h = legend('$\bar{\theta}$','$\theta(p_0)$','$\theta(p_{opt})$');
set(h,'Interpreter','latex','FontSize',fsz);
xlabel('time', 'FontSize',fsz);
ylabel('val [rad]', 'FontSize',fsz);
title('$\theta$ for different parameters','interpreter','latex','FontSize',fsz);


end
