function [] = bobEx3p2 ()

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
courant = [0.8, 0.4, 0.2, 0.1];
for i = 1:length(courant)
    rms(i) = solver(init_fun, courant(i), phase_speed, dx, end_time);
end

display((log(rms(end-1)) - log(rms(end))) / (log(courant(end-1)) - log(courant(end))))

figure;
loglog(courant(2:end), rms(2:end), 'linewidth', 2.0,'color','k');
xlabel('Courant number','fontsize',15)
ylabel('RMSE','fontsize',15)
title('dx = 1/20','fontsize',15)
set(gca,'fontsize',15)

% Part c)

dx = 1/80;
courant = [0.8, 0.4, 0.2, 0.1];
for i = 1:length(courant)
    rms(i) = solver(init_fun, courant(i), phase_speed, dx, end_time)
end

display((log(rms(end-1)) - log(rms(end))) / (log(courant(end-1)) - log(courant(end))))

figure;
loglog(courant, rms, 'linewidth', 2.0,'color','k');
xlabel('Courant number','fontsize',15)
ylabel('RMSE','fontsize',15)
title('dx = 1/80','fontsize',15)
set(gca,'fontsize',15)


end

function rms = solver(init_fun, courant, phase_speed, dx, end_time)

x = (0:dx:1-dx)';
dt = courant * dx / phase_speed;

LELE_MAT = construct_lele_mat(length(x));

ASSELIN_PARAM = 0.1;

state_filtered = init_fun(x);

num_steps = max(round(end_time/dt), 2);
step_num = 0;

% Take one midpoint step to predict current state.
state_midpoint = state_filtered + 0.5*dt*compute_tendency(state_filtered, phase_speed, dx, LELE_MAT);
state_current = state_filtered + dt*compute_tendency(state_midpoint, phase_speed, dx, LELE_MAT);

% Predict next state using leapfrog.
state_next = state_filtered + 2*dt*compute_tendency(state_current, phase_speed, dx, LELE_MAT);
step_num = step_num + 1;

% Leapfrog
while step_num < num_steps

    state_filtered = state_current + ASSELIN_PARAM * (state_next - 2*state_current + state_filtered);
    state_current = state_next;
    state_next = state_filtered + 2*dt*compute_tendency(state_next, phase_speed, dx, LELE_MAT);
    
    step_num = step_num + 1;
    
end

exact = init_fun(x - phase_speed*end_time);
rms = sqrt(mean((state_current - exact).^2));

end

function LELE_MAT = construct_lele_mat(num_points)

LELE_MAT = spalloc(num_points, num_points, 3*num_points);

LELE_MAT(1,1) = 14;
LELE_MAT(1,2) = 5;
LELE_MAT(1,end) = 5;

for i = 2:num_points-1
   LELE_MAT(i,i-1) = 5;
   LELE_MAT(i,i+1) = 5;
   LELE_MAT(i,i) = 14;
end

LELE_MAT(end,end) = 14;
LELE_MAT(end,1) = 5;
LELE_MAT(end,end-1) = 5;

LELE_MAT = LELE_MAT / 24;

end

function dt_solution = compute_tendency(solution, phase_speed, dx, LELE_MAT)

F = zeros(length(solution) + 4, 1);
F(1:2) = solution(end-1:end);
F(end-1:end) = solution(1:2);
F(3:end-2) = solution;

lele_rhs = (11*(F(4:end-1) - F(2:end-3))/(2*dx) + (F(5:end) - F(1:end-4))/(4*dx)) / 12;

dx_solution = LELE_MAT \ lele_rhs;

dt_solution = -phase_speed * dx_solution;


end