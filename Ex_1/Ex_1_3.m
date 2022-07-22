clear all; close all; clc;

%% Lotka Volterra
% Optimization problem

load HL_population.mat
global HL
global H
global L
global time

H = H/1e3; L = L/1e3;

p0 = [1.1; 0.4; 0.3; 0.4; 20; 32];

fn = @lv_loss;
lb = [0 0 0 0 7 9];
ub = [3 1 4 1 150 150];
[p_opt,res] = fmincon(fn,p0,[],[],[],[],lb,ub)


disp(res)
disp('OPTIMAL PARAMETERS:')
disp({'b', 'p', 'r', 'd'})
disp(p_opt')


% Final solution optimized
options = odeset('AbsTol', 10^-12, 'RelTol', 10^-12, 'MaxStep', 0.1);
lv_opt = @(t,u) lvfun(p_opt(1),p_opt(2),p_opt(3),p_opt(4),t,u);     

[t_fin,x_fin] = ode45(lv_opt , [0:0.1:30] , [p_opt(5); p_opt(6)]);
year_f = years(1) + t_fin*2;



% figure(1)

plot(years,H.*1e3,'r--o','LineWidth',.5)

hold on
plot(years,L.*1e3,'b--o','LineWidth',.5)
plot(year_f,x_fin(:,1).*1e3,'r','Linewidth',2)
plot(year_f,x_fin(:,2).*1e3,'b','Linewidth',2)

ylabel ('Population State')
xlabel ('Years')
title ('Lotka-Volterra Model')
legend ('Hare','Lynx','LV Hare','LV Lynx');
xlim([1845 1903])

grid on

function loss = lv_loss(par)

global HL
global H
global L
global time

b = par(1);
p = par(2);
r = par(3);
d = par(4);

H0 = par(5);
L0 = par(6);

timespan = [1:0.1:30];

f = @(t,u) lvfun(b, p, r, d, t, u);
[t_mod,x_mod] = ode45(f,timespan,[H0;L0]);

H_mod = x_mod(:,1);
L_mod = x_mod(:,2);

H_data = interp1(time,H,t_mod);
L_data = interp1(time,L,t_mod);

loss = sqrt(sum((H_data-H_mod).^2 + (L_data-L_mod).^2));

end

