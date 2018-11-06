function [] = bobEx3p1 ()

% Part a)

init_fun = @(x) sin(2*pi*x).^6;
phase_speed = 0.1;
end_time = 50;
courant = 0.1;

dx = [1/20, 1/40, 1/80, 1/160];
for i = 1:length(dx)
    rms(i) = solver(init_fun, courant, phase_speed, dx(i), end_time);
end

display((log(rms(end-1)) - log(rms(end))) / (log(dx(end-1)) - log(dx(end))))
fourth_order_line = exp(4*(log(dx) - log(dx(end))) + log(rms(end)));

figure;
loglog(dx, fourth_order_line, 'linewidth', 2.0, 'linestyle', '--','color','k');
hold all
loglog(dx, rms, 'linewidth', 2.0,'color','k');
xlabel('dx','fontsize',15)
ylabel('RMSE','fontsize',15)
set(gca,'fontsize',15)
legend('O(dx^4)', 'Actual','location','northwest')

% Part b)

dx = 1/20;
courant = [0.1, 0.2, 0.4, 0.8];
for i = 1:length(courant)
    rms(i) = solver(init_fun, courant(i), phase_speed, dx, end_time);
end

figure;
loglog(courant(1:3), rms(1:3), 'linewidth', 2.0,'color','k');
xlabel('Courant number','fontsize',15)
ylabel('RMSE','fontsize',15)
title('dx = 1/20','fontsize',15)
set(gca,'fontsize',15)

% Part c)

dx = 1/80;
courant = [0.1, 0.2, 0.4, 0.8];
for i = 1:length(courant)
    rms(i) = solver(init_fun, courant(i), phase_speed, dx, end_time);
end

figure;
loglog(courant(1:3), rms(1:3), 'linewidth', 2.0,'color','k');
xlabel('Courant number','fontsize',15)
ylabel('RMSE','fontsize',15)
title('dx = 1/80','fontsize',15)
set(gca,'fontsize',15)


end

function rms = solver(init_fun, courant, phase_speed, dx, end_time)

x = (0:dx:1-dx)';
dt = courant * dx / phase_speed;

solution = init_fun(x);
tendency = zeros(length(x),3);
CURRENT = 3; PREVIOUS = 2; TWO_STEPS_AGO = 1;

num_steps = max(round(end_time/dt), 2);
step_num = 0;

% Take two RK3 steps
for i = 1:2
    tendency1 = compute_tendency(solution, phase_speed, dx);
    tendency(:,i) = tendency1; % Tendencies required by the Adams-Bashforth 
    
    pred_sol1 = solution + dt * tendency1 / 3;
    tendency2 = compute_tendency(pred_sol1, phase_speed, dx);
    
    pred_sol2 = pred_sol1 + dt * (tendency2 - (5/9)*tendency1) * (15/16);
    tendency3 = compute_tendency(pred_sol2, phase_speed, dx);
    
    solution = pred_sol2 + dt * (tendency3 - (153/128)*(tendency2 - (5/9)*tendency1)) * (8/15);
    
    step_num = step_num + 1;
end

% AB3
while step_num < num_steps
    tendency(:,CURRENT) = compute_tendency(solution, phase_speed, dx);
    
    solution = solution + dt * (23*tendency(:,CURRENT) - 16*tendency(:,PREVIOUS) + 5*tendency(:,TWO_STEPS_AGO)) / 12;
    
    tendency(:,TWO_STEPS_AGO) = tendency(:,PREVIOUS);
    tendency(:,PREVIOUS) = tendency(:,CURRENT);
        
    step_num = step_num + 1;
    
end

exact = init_fun(x - phase_speed*end_time);
rms = sqrt(mean((solution - exact).^2));

end

function dt_solution = compute_tendency(solution, phase_speed, dx)

F = zeros(length(solution) + 4, 1);
F(1:2) = solution(end-1:end);
F(end-1:end) = solution(1:2);
F(3:end-2) = solution;

dt_solution = -phase_speed * (F(1:end-4) - 8*F(2:end-3) + 8*F(4:end-1) - F(5:end)) / (12*dx);


end