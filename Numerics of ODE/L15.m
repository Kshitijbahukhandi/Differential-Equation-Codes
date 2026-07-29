%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Forward Euler Method on a Variable Mesh
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
%     Forward Euler Method with Variable Mesh
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
  ETA=1.2;
  disp('Variable Mesh Numerical Solution for the system of IVPs Using Forward Euler Method');
  disp('    ');
  fprintf('N=%4.2f\n',N);
  disp('  ');
  for I=1:8
    % Formation of variable mesh lengths
       SUM=1.0;
       for L=1:N-1
           SUM=SUM+ETA^L;
       end
       H(2)=1.0/SUM;
       for L=2:N
           H(L+1)=ETA*H(L);
       end 
  % Formation of grid points
    T(1)=0.0;
    for L=2:N+1
          T(L)=T(L-1)+H(L);   
    end 
    T(N+1)=1.0;
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
           H0=H(L+1);
       U(L+1)=U(L)+H0*(3*U(L)-V(L)); 
       V(L+1)=V(L)+2*H0*U(L);
       end 
  % Calculate Errors
        ERU=abs(U(N+1)-EU(N+1));
        ERV=abs(V(N+1)-EV(N+1));
     % Print errors
       fprintf('N=%4.2f, ERU=%10.4e, ERV=%10.4e\n',N,ERU, ERV);
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
    figure (2);
  plot(T,V,'-g',T,EV,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values v');
  title ('Exact vs Numerical Solution')
  
