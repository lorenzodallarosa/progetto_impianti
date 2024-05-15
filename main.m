v_s   = 3;      % m/s
M     = 1e4;    % kg
g     = 9.81;   % m/s^2
A_h   = 6e-3;   % m^2
c_pn  = 8e3;    % Ns/m
k_pn  = 5e6;    % N/m
h_0   = 0.75;   % m
A_g   = 4.8e-3; % m^2
P_g0  = 4e6;    % Pa
n     = 1.3;    % pure number
rho   = 850;    % kg/m^3
k_s   = 2;      % pure number
A_s0  = 1e-4;   % m^2
c_max = 0.5;    % m
B_e   = 800e6;  % Pa
F_0   = 25e3;   % N
tau   = 3;      % s

%% Variables Vector

% u = [z; w; s; P_h] dove w = z';

%% Functions

function force = F_up(time)
force = F_0 * exp(- time / tau);
end

function pressure = P_g(u)
if u(4) >= P_g0
  pressure = P_g0 * (1 - (A_h) / (A_g * h_0) * (u(1) - u(3))).^(-n);
else
  pressure = P_g0;
end
end

function area = A_s(u)
area = A_s0 * (1 - ((u(1) - u(3)) / c_max).^2);
end

function flow_rate = Q_s(u)
if u(4) >= P_g0
  flow_rate = A_s(u) * sqrt((2 * abs(u(4) - P_g(u))) / (rho * k_s)) * sign(u(4) - P_g(u));
else
  flow_rate = 0;
end
end

function u_d = f(u, time)
u_d = ones(4, 1);
u_d(1) = u(2);
u_d(2) = g - (A_h / M) * u(4) - F_up(time) / M;
u_d(3) = (A_h / c_pn) * u(4) - (k_pn / c_pn) * u(3);
u_d(4) = B_e * (A_h * (u_d(1) - u_d(3)) - Q_s(u)) / (A_h * (c_max - u(1) + u(3)));
end

