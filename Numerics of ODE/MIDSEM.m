%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a Second-Order Initial Value Problem
% using the Fourth-Order Runge-Kutta Method
%
% Differential Equation:
%     u'' = 6u^2
%
% Problem Type:
%     Second-Order Initial Value Problem (IVP)
%
% Numerical Method:
%     Classical Fourth-Order Runge-Kutta (RK4) Method
%
% Language:
%     MATLAB
%
% Developed as part of the coursework for
% Numerical Methods for Ordinary Differential Equations
% M.Sc. Applied Mathematics
% South Asian University
%
% Author: Kshitij Bahukhandi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = @(t, u, v) 6*u^2;

% Define the initial conditions
t0 = 0;
u0 = 1;
v0 = -2;

% Define the exact solution
exact_solution = @(t) 1 / (1 + t)^2;

% Define the mesh sizes
h_values = [1/40, 1/80, 1/160, 1/320];

% Define the points where we want to compute absolute errors
r_values = [0.5, 1.0];

% Initialize arrays to store absolute errors
errors = zeros(length(h_values), length(r_values));

% Perform RK4 method for each mesh size
for i = 1:length(h_values)
    h = h_values(i);
    for j = 1:length(r_values)
        r = r_values(j);
        
        % Initialize arrays for RK4 method
        t = t0:h:1;
        N = length(t);
        U = zeros(1, N);
        V = zeros(1, N);
        
        % Set initial values
        U(1) = u0;
        V(1) = v0;
        
        % Perform RK4 iterations
        for k = 1:N-1
            k1 = h * V(k);
            l1 = h * f(t(k), U(k), V(k));
            
            k2 = h * (V(k) + 0.5*l1);
            l2 = h * f(t(k) + 0.5*h, U(k) + 0.5*k1, V(k) + 0.5*l1);
            
            k3 = h * (V(k) + 0.5*l2);
            l3 = h * f(t(k) + 0.5*h, U(k) + 0.5*k2, V(k) + 0.5*l2);
            
            k4 = h * (V(k) + l3);
            l4 = h * f(t(k) + h, U(k) + k3, V(k) + l3);
            
            U(k+1) = U(k) + (k1 + 2*k2 + 2*k3 + k4) / 6;
            V(k+1) = V(k) + (l1 + 2*l2 + 2*l3 + l4) / 6;
        end
        
        % Compute absolute error at the specified point
        exact_value = exact_solution(r);
        approx_value = U(round(r / h) + 1); % Approximation at r
        errors(i, j) = abs(exact_value - approx_value);
    end
end

% Display the absolute errors for each mesh size and point
disp('Absolute Errors:');
disp(errors);
% This code defines the differential equation, initial conditions, exact solution, mesh sizes, and points where absolute errors are computed. It then performs the RK4 method for each mesh size and computes the absolute errors. The results are stored in the errors array and displayed at the end.