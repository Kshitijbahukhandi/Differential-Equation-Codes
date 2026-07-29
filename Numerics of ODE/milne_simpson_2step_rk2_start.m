%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Milne-Simpson Method
%
% Differential Equation:
%     du/dt = -2*t*u^2
%
% Exact Solution:
%     u(t) = 1/(1+t^2)
%
% Initialization:
%     Second-Order Classical Explicit Runge-Kutta (RK2) Method
%
% Numerical Method:
%     Two-Step Milne-Simpson Method
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

clear all; 
for I = 1:8 
    N = 10; 
    H = 1.0 / N; 
    disp('Numerical Solution for IVP using Milne-Simpson Method'); 
    disp('  '); 
    % Formation of Nodal Points  
    for L = 1:N+1 
        T(L) = (L-1) * H; 
    end 
     
    % Exact Solution & RHS Function  
    for L = 1:N+1 
        EU(L) = 1.0 / (1 + T(L) * T(L)); 
    end 
     
    % Initial Conditions 
    U(1) = 1.0; 
    tic % CPU time Starts 
    % Solution at t=h by the help of second order classical explicit RK-Method 
    K1 = -2 * H * T(1) * U(1) * U(1); 
    K2 = -2 * H * (T(1) + H) * (U(1) + K1) * (U(1) + K1); 
    U(2) = U(1) + (K1 + K2) / 2; 
    % Solution at advanced Time 
    for L = 2:N 
        U0 = U(L); 
        UU0 = U(L-1); 
        T0 = T(L); 
        TT0 = T(L-1); 
        U(L+1) = UU0 - 4 * H * T0 * U0 * U0; 
    end 
    % Calculate Errors 
    ERU = abs(U(N+1) - EU(N+1)); 
     
    % Print Errors 
    fprintf('N=%4.2f, ERU=%10.4e\n', N, ERU); 
    disp(' '); 
    N = 2 * N; 
end

toc % CPU Time ends 
% Elapsed time is 0.019099 seconds. 
% Plot the Results 
figure(1); 
plot(T, U, '-g', T, EU, '--r'); 
xlabel('t-values'); 
ylabel('Numerical and Exact values u'); 
title('Exact vs Numerical Solution');

