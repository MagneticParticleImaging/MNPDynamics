clear
close all
tic

addpath('../matlab')
load('../matlab/FV_meshes/mesh_4.mat')

params.kAnis = 11000;
%params.n = [0;0;1];
f = 25000;
%params.n = @(t) [cos(2*pi*f*t);sin(2*pi*f*t);0];
B =@(t)0.012*[cos(.34*2*pi*f*t)-.05;sin(.7*2*pi*f*t)+.1;sin(2*pi*f*t)];
params.p1 = 0;
params.p3 = 0;
%B = @(t)0.012*[0;0;sin(2*pi*f*t)];

t = linspace(0,4/f, 1000);

[t, exp, y] = simulation_FV(B, t, tr, params);

xexp = exp(:,1);
yexp = exp(:,2);
zexp = exp(:,3);

figure
subplot(1,2,1)
plot(t, [xexp,yexp,zexp])
legend({'$\bar{m}_x$','$\bar{m}_y$','$\bar{m}_z$'}, 'Interpreter','latex')
subplot(1,2,2)
dt = diff(t);
dxdt = diff(xexp)./dt;
dydt = diff(yexp)./dt;
dzdt = diff(zexp)./dt;
plot(t(1:end-1),[dxdt, dydt, dzdt]) 
legend({'$\frac{\partial \bar{m}_x}{\partial t}$','$\frac{\partial \bar{m}_y}{\partial t}$','$\frac{\partial \bar{m}_z}{\partial t}$'}, 'Interpreter','latex')
% Note: the result has to be multiplied by M_S*V_C to obtain the magnetic
% moment

% plot the probability distribution over time

% figure
% for i=1:length(t)
%     trisurf(tr.fMat, tr.vMat(:,1), tr.vMat(:,2), tr.vMat(:,3), y(i,:), 'EdgeColor', 'none')
%     title(num2str(i/length(t)))
%     caxis([min(min(y)), max(max(y))]);
%     colorbar()
%     drawnow()
% end
% 
