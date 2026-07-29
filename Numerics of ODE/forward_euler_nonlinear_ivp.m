%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a Nonlinear Initial Value Problem
% using the Forward Euler Method
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
%     Forward Euler Method
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
  disp('Numerical Solution for IVP Using Forward Euler Method');
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
          U(L+1)=U(L)-2*H*T(L)*U(L)*U(L);
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
  figure(4);
  plot(T,U,'b+'); hold on
  plot(T,EU,'-r'); hold off
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title('Numerical vs Exact Solution')

