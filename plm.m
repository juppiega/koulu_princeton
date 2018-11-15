function [] = plm (dx, USE_STRANG_SPLITTING, USE_EASTER)
% plm (dx, USE_STRANG_SPLITTING, USE_EASTER)
%
% Arguments:
% dx                    : (OPTIONAL) Grid spacing (Default: 0.02)
% USE_STRANG_SPLITTING  : (OPTIONAL) Whether to use Strang splitting (Default: true)
% USE_EASTER            : (OPTIONAL) Whether to use Easter's pseudocompressibility approach (Default: true)

if nargin == 0
    dx = 0.02;
end
if nargin <= 1
    USE_STRANG_SPLITTING = true;
end
if nargin <= 2
    USE_EASTER = true;
end

% Initialization
[X, Y, tracer, dt] = initialize(dx);
plot_tracer(X, Y, tracer, 0)

% Advect 2.5 units and plot
tracer = advect_tracer(X, Y, tracer, dx, dt, 0, 2.5, USE_STRANG_SPLITTING, USE_EASTER);
plot_tracer(X, Y, tracer, 2.5)

% Advect from time 2.5 to 5 and plot
tracer = advect_tracer(X, Y, tracer, dx, dt, 2.5, 5, USE_STRANG_SPLITTING, USE_EASTER);
plot_tracer(X, Y, tracer, 5)


end

% ===============================================================================================
function tracer = advect_tracer(X, Y, tracer, dx, dt, time_begin, time_end, USE_STRANG_SPLITTING, USE_EASTER)

time = time_begin;
first_x_then_y = true;

% ******************************************************************************
% MAIN LOOP
% ******************************************************************************
while time < time_end

    % Adapt timestep at the very last step, if time_end is not a multiple of dt
    if time_end - time < dt
        dt = time_end - time;
    end
    
    % Compute interface winds, centered in time.
    [u, v] = described_wind(time + dt/2, X, Y, dx);
    
    rho = 1; % Pseudodensity
    tracer_rho = tracer(2:end-1, 2:end-1) .* rho; % Tracer times density
    
    if ~USE_STRANG_SPLITTING
        first_x_then_y = true;
    end
    
    % Advance in time by reversing the order of integration every timestep (Strang splitting)
    if first_x_then_y
        [tracer, rho, tracer_rho] = compute_x_direction(dx, dt, u, tracer, rho, tracer_rho, USE_EASTER);    
        tracer                    = compute_y_direction(dx, dt, v, tracer, rho, tracer_rho, USE_EASTER);
        first_x_then_y = false;
    else
        [tracer, rho, tracer_rho] = compute_y_direction(dx, dt, v, tracer, rho, tracer_rho, USE_EASTER);
        tracer                    = compute_x_direction(dx, dt, u, tracer, rho, tracer_rho, USE_EASTER); 
        first_x_then_y = true;
    end
    
    time = time + dt;
    
end

end

% ==========================================================================================
function [tracer, rho, tracer_rho] = compute_x_direction(dx, dt, u, tracer, rho, tracer_rho, USE_EASTER)

% Compute one-dimensional scalar flux
F_x = compute_fluxes(dx, dt, u, tracer);

if USE_EASTER
    % Update pseudodensity
    rho = rho - (dt/dx) * (u(:,2:end) - u(:,1:end-1)); % Durran's equation (5.95)
    % Update tracer mass
    tracer_rho = tracer_rho - (dt/dx) * (F_x(:,2:end) - F_x(:,1:end-1)); % eq. (5.96)
    % Update tracer concentrations
    tracer(2:end-1, 2:end-1) = tracer_rho ./ rho; % eq. (5.97)
else
    tracer(2:end-1, 2:end-1) = tracer(2:end-1, 2:end-1) - (dt/dx) * (F_x(:,2:end) - F_x(:,1:end-1));
end

end

% ==========================================================================================
function [tracer, rho, tracer_rho] = compute_y_direction(dx, dt, v, tracer, rho, tracer_rho, USE_EASTER)

% With some transpose-magic, the same flux-computation function can be used for both u and v
F_y = compute_fluxes(dx, dt, v', tracer')';

if USE_EASTER
    rho = rho - (dt/dx) * (v(2:end,:) - v(1:end-1,:)); % Durran's equation (5.98)
    tracer_rho = tracer_rho - (dt/dx) * (F_y(2:end,:) - F_y(1:end-1,:)); % eq. (5.99)
    tracer(2:end-1, 2:end-1) = tracer_rho ./ rho; % eq. (5.100)
else
    tracer(2:end-1, 2:end-1) = tracer(2:end-1, 2:end-1) - (dt/dx) * (F_y(2:end,:) - F_y(1:end-1,:));
end

end

% ============================================
function slope = compute_limited_slope(tracer)
% This function selects the smallest magnitude slope from a set of three slopes 
% (first order forward and backward, and second order centered).

N = size(tracer,1);
slope_centered = zeros(N-2,N-2);
slope_right = zeros(N-2,N-2);
slope_left = zeros(N-2,N-2);

% Compute the three slopes
slope_centered = (tracer(2:end-1,3:end) - tracer(2:end-1,1:end-2)) / 2;
slope_right = 2 * (tracer(2:end-1,3:end) - tracer(2:end-1,2:end-1));
slope_left = 2 * (tracer(2:end-1,2:end-1) - tracer(2:end-1,1:end-2));

slope_interior_points = zeros(N-2,N-2);
p_i = slope_centered > 0 & slope_left > 0 & slope_right > 0; % All slopes positive
n_i = slope_centered < 0 & slope_left < 0 & slope_right < 0; % All negative
slope_interior_points(p_i) = min([slope_centered(p_i), slope_left(p_i), slope_right(p_i)], [], 2);
slope_interior_points(n_i) = max([slope_centered(n_i), slope_left(n_i), slope_right(n_i)], [], 2);
% Now, local minima and maxima have slope = 0, corresponding to an upwind scheme.

slope = zeros(size(tracer));
slope(2:end-1,2:end-1) = slope_interior_points;

end

% =================================================
function F = compute_fluxes(dx, dt, u_or_v, tracer)
% Computes one-dimensional fluxes in either direction.

c = u_or_v * dt/dx; % Courant

slope = compute_limited_slope(tracer);

% Use different formulae for positive and negative flow directions. 
% The first one is from the lecture notes, and the second one is straightforward to derive.
F_positive_c = u_or_v .* (tracer(2:end-1,1:end-1) - 0.5 * slope(2:end-1,1:end-1) .* (c - 1));
F_negative_c = u_or_v .* (tracer(2:end-1,2:end) - 0.5 * slope(2:end-1,2:end) .* (c + 1));
F = zeros(size(u_or_v));
F(c <  0) = F_negative_c(c <  0);
F(c >= 0) = F_positive_c(c >= 0);

end

% ===========================================
function [u, v] = described_wind(t, X, Y, dx)
% Described wind components

u = zeros(size(X,1)-2, size(X,2)-1);
v = zeros(size(X,1)-1, size(X,2)-2);

% Sample at interfaces.
X_i = X + 0.5*dx;
Y_i = Y + 0.5*dx;

% Described velocities
u = sin(pi * X_i(2:end-1, 1:end-1)).^2 .* sin(2*pi * Y_i(2:end-1, 1:end-1)) * cos(pi*t/5);
v = -sin(2*pi * X_i(1:end-1, 2:end-1)) .* sin(pi * Y_i(1:end-1, 2:end-1)).^2 * cos(pi*t/5);

end

% ===========================================
function [X, Y, tracer, dt] = initialize(dx)
% Initialization function

MAX_COURANT = 0.5;

x = 0:dx:1;
y = x;
[X,Y] = meshgrid(x, y); % Create 2D meshes of coordinates for convenience.

% Gaussian bell
r = min(1, 4*sqrt((X-1/4).^2 + (Y-1/4).^2));
tracer = 0.5 * (1 + cos(pi*r));

% Compute timestep based on maximum L1-norm of the wind.
[u0, v0] = described_wind(0, X, Y, dx);
U_L1_norm = abs(u0(:)) + abs(v0(:));
dt = MAX_COURANT * dx / max(U_L1_norm);

end

% ===================================
function plot_tracer(X, Y, tracer, t)

% Plot
figure;
surf(X, Y, tracer, 'edgecolor','none')
view(2)
title(['Tracer concentration @ t = ', num2str(t)], 'fontsize', 15)
set(gca, 'fontsize', 15)
set (gca, 'color', 'b')
ylabel('Y', 'fontsize', 15)
xlabel('X', 'fontsize', 15)
colorbar('fontsize', 15)
caxis([0,1])

% Analyze
format compact
display(['Maximum value @ t = ', num2str(t), ': ', num2str(max(tracer(:)))])
cell_area = (X(1,2) - X(1,1)) * (Y(2,1) - Y(1,1));
total_tracer_amount = cell_area * sum(tracer(:));
display(['Total amount @ t = ', num2str(t), ': ', num2str(total_tracer_amount)])
fprintf('\n')

end
