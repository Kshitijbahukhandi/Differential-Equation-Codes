%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Backward Euler Method
%
% System of Differential Equations:
%     du/dt = 3u - v
%     dv/dt = 2u
%
% Exact Solutions:
%     u(t) = exp(t)
%     v(t) = 2*exp(t)
%
% Numerical Method:
%     Backward Euler Method
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%cal Solution of system of IVPs - Backward Euler Method

  clear all;  % clear all variables in memory
  N=10; 
  for I=1:8
  H=1.0/N;
  disp('Numerical Solution of system of IVPs Using Backward Euler Method');
  disp('    ');
  fprintf('N=%4.2f, H=%6.4f\n',N,H);
  disp('  ');
  % Formation of Nodal Points
    for L=1:N+1
          T(L)=(L-1)*H;   
    end 
    % Exact Solution & RHS function
    for L=1:N+1
        EU(L)= exp(T(L));
        EV(L)=2*exp(T(L));
    end
    % Initial conditions
 %   U(1)=EU(1);
 %   V(1)=EV(1);
    U(1)=1.0;
    V(1)=2.0;
     tic  % CPU time starts
    % Solution at advanced time
       for L=1:N   
     DET=1-3*H+2*H*H;
     DETU=U(L)-H*V(L);
     DETV=(1-3*H)*V(L)+2*H*U(L);
     U(L+1)=DETU/DET;
     V(L+1)=DETV/DET;
       end 
  % Calculate Errors
        ERU=abs(U(N+1)-EU(N+1));
        ERV=abs(V(N+1)-EV(N+1));
     % Print errors
       fprintf('N=%4.2f, ERU=%10.4e, ERV=%10.4e\n',N,ERU,ERV);
  disp('  ');
  N=2*N;
  end
   toc  % CPU time ends
    % Plot the results
  figure (1);
  plot(T,U,'-g',T,EU,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title ('Exact vs Numerical Solution')
  figure(2);
  plot(T,V,'-g',T,EV,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values v');
  title ('Exact vs Numerical Solution')

