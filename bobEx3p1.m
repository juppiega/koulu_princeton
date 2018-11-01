function [state_history] = lorenz (X0, r, dt, end_time)
% [state_history] = lorenz (X0, r, dt, end_time)
%
% Arguments:
% X0 : Initial position, e.g. X0 = [10, 10, 10]
% r  : Constant parameter, determining colution characteristics (~25 for chaos)
% dt : Timestep (0.03 recommended)
% end_time : Length of the simulation in units of time
%
% Output:
% state_history : Full history of the model state (dimension: 3 x [NUM_STEPS] array)


if isrow(X0); X0 = X0; end

state_history = solver(X0, r, dt, end_time);

figure;
plot(state_history(1,:), state_history(2,:),'.'); 
xlabel('x')
ylabel('y')

figure;
plot(state_history(1,:), state_history(3,:),'.');
xlabel('x')
ylabel('z')

end

function state_history = solver(X0, r, dt, end_time)

X = X0;
tendency = zeros(3,3);
CURRENT = 3; PREVIOUS = 2; TWO_STEPS_AGO = 1;
num_steps = max(ceil(end_time/dt), 2);
state_history = zeros(3, num_steps); state_history(:,1) = X;
step_num = 1;

% Take two RK3 steps
for i = 1:2
    tendency1 = computeDX(X, r);
    tendency(:,i) = tendency1; % Tendencies required by the Adams-Bashforth 
    
    pred_X1 = X + dt * tendency1 / 3;
    tendency2 = computeDX(pred_X1, r);
    
    pred_X2 = pred_X1 + dt * (tendency2 - (5/9)*tendency1) * (15/16);
    tendency3 = computeDX(pred_X2, r);
    
    X = pred_X2 + dt * (tendency3 - (153/128)*(tendency2 - (5/9)*tendency1)) * (8/15);
    
    step_num = step_num + 1;
    state_history(:, step_num) = X;
end

% AB3
while step_num < num_steps
    tendency(:,CURRENT) = computeDX(X, r);
    
    X = X + dt * (23*tendency(:,CURRENT) - 16*tendency(:,PREVIOUS) + 5*tendency(:,TWO_STEPS_AGO)) / 12;
    
    tendency(:,TWO_STEPS_AGO) = tendency(:,PREVIOUS);
    tendency(:,PREVIOUS) = tendency(:,CURRENT);
    
    step_num = step_num + 1;
    state_history(:, step_num) = X;
end

end

function dX = computeDX(X, r)

dX = zeros(3,1);
dX(1) = -3*(X(1)-X(2));
dX(2) = -X(1)*X(3) + r*X(1) - X(2);
dX(3) = X(1)*X(2) - X(3);

end