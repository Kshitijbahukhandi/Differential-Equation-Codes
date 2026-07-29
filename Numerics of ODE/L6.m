%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Exponential Approach
%
% Differential Equation:
%     du/dt = lambda*u
%
% Exact Solution:
%     u(t) = exp(lambda*t)
%
% Numerical Method:
%     Exponential Approach Method
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
  LAMDA=-10.0;
  N=10; 
  for I=1:8
  H=1.0/N;
  disp('Numerical Solution for IVP Using Exponential Approach');
  disp('    ');
  fprintf('N=%4.2f, H=%6.4f\n',N,H);
  disp('  ');
  % Formation of Nodal Points
    for L=1:N+1
          T(L)=(L-1)*H;   
    end 
    % Exact Solution & RHS function
    for L=1:N+1
        EU(L)= exp(LAMDA*T(L)); 
    end
    % Initial conditions
 %   U(1)=EU(1);
    U(1)=1.0;
     tic  % CPU time starts
    % Solution at advanced time
       for L=1:N
           C=LAMDA*H;
           % First Order Approx
             D=1+C;
             % Second Order Approx
             D=1+C+0.5*C*C;
             % Third Order Approx
             D=1+C+0.5*C*C+(C*C*C/6);
             % Fourth Order Approx
             D=1+C+0.5*C*C+(C*C*C/6)+(C*C*C*C/24);
    %----
       U(L+1)=D*U(L);
       end 
  % Calculate Errors
        ERU=abs(U(N+1)-EU(N+1));
     % Print errors
       fprintf('N=%4.2f, ERU=%10.4e\n',N,ERU);
  disp('  ');
  N=2*N;
  end
   toc  % CPU time ends
    % Plot t.he results
  figure (3);
  plot(T,U,'--b',T,EU,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title ('Exact vs Numerical Solution')
  figure(4);
  plot(T,U,'b+'); hold on
  plot(T,EU,'-r'); hold off
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title('Numerical vs Exact Solution')
