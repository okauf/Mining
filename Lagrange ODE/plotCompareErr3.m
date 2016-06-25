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
else
    load 'data/compareErr3.mat';
end

% object names:
% absErr, relErr, absSTErr, relSTErr, fctVal, p_0, maxsdiff, maxtdiff

% font size
fsz = 15;

figure;
%subplot(2,2,1);
semilogy(squeeze(fctVal(1,1,:)),'-k');
hold on;
semilogy(squeeze(fctVal(2,1,:)),'-.k');
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error', 'FontSize',fsz);
title('$0.5 \| \Phi (\bar{x},u,p) \|^2$','interpreter','latex', 'FontSize',fsz);

figure;
%subplot(2,2,2);
semilogy(squeeze(absErr(1,1,:)),'-k');
hold on;
semilogy(squeeze(absErr(2,1,:)),'-.k');
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error', 'FontSize',fsz);
title('$0.5 \| \bar{x}-x(p)\|^2$','interpreter','latex', 'FontSize',fsz);

figure;
%subplot(2,2,3);
semilogy(squeeze(maxsdiff(1,1,:)),'-k');
hold on;
semilogy(squeeze(maxsdiff(2,1,:)),'-.k');
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error [m]', 'FontSize',fsz);
title('$\max(|\bar{s}_n - s_n(p)|)$','interpreter','latex', 'FontSize',fsz);

figure;
%subplot(2,2,4);
semilogy(squeeze(maxtdiff(1,1,:)),'-k');
hold on;
semilogy(squeeze(maxtdiff(2,1,:)),'-.k');
legend('interior-point','sqp');
xlabel('iter', 'FontSize',fsz);
ylabel('error [rad]', 'FontSize',fsz);
title('$\max(|\bar{\theta}_n - \theta_n(p)|)$','interpreter','latex', 'FontSize',fsz);

end
