%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a Nonlinear Initial Value Problem
% using the Second-Order Explicit Runge-Kutta Method
%
% Differential Equation:
%     du/dt = -2*t*u^2
%
% Exact Solution:
%     u(t) = 1/(1+t^2)
%
% Problem Type:
%     Nonlinear Initial Value Problem (IVP)
%
% Numerical Method:
%     Second-Order Explicit Runge-Kutta (RK2) Method
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

clear all;  % clear all variables in memory
  N=10; 
  for I=1:8
  H=1.0/N;
  disp('Numerical Solution for IVP second order explicit RK Method');
  disp('    ');
  fprintf('N=%4.2f, H=%6.4f\n',N,H);
  disp('  ');
  % Formation of Nodal Points
    for L=1:N+1
          T(L)=(L-1)*H;   
    end 
    % Exact Solution & RHS function
    for L=1:N+1
        EU(L)= 1.0/(1+T(L)*T(L)); 
    end
    % Initial conditions
 %   U(1)=EU(1);
    U(1)=1.0;
     tic  % CPU time starts
    % Solution at advanced time
       for L=1:N 
           T0=T(L); U0=U(L);
       %---Classical (Heun) Method
           K1= -2*H*T0*U0*U0;
           K2= -2*H*(T0+H)*(U0+K1)*(U0+K1);
           U(L+1)=U(L)+(K1+K2)/2;
       %---Optimal Method
           K1= -2*H*T0*U0*U0;
           K2= -2*H*(T0+2*H/3)*(U0+2*K1/3)*(U0+2*K1/3);
           U(L+1)=U(L)+(K1+3*K2)/4;
       %---Euler Tangent Method 
           K1= -2*H*T0*U0*U0;
           K2= -2*H*(T0+H/2)*(U0+K1/2)*(U0+K1/2);
           U(L+1)=U(L)+K2;
       end 
  % Calculate Errors
        ERU=abs(U(N+1)-EU(N+1));
     % Print errors
       fprintf('N=%4.2f, ERU=%10.4e\n',N,ERU);
  disp('  ');
  N=2*N;
  end
   toc  % CPU time ends
    % Plot the results
  figure(1);
  plot(T,U,'--g'); 
  xlabel('t-values');
  ylabel('Numerical values u');
  title('Numerical Solution')
  figure (2);
  plot(T,EU,'--r');
  xlabel('t-values');
  ylabel('Exact values u');
  title ('Exact Solution')
  figure (3);
  plot(T,U,'-g',T,EU,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title ('Exact vs Numerical Solution')
