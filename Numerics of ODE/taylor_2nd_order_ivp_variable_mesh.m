%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Second-Order Taylor Series Method
% on a Variable Mesh
%
% Differential Equation:
%     du/dt = -2*t*u^2
%
% Exact Solution:
%     u(t) = 1/(1+t^2)
%
% Numerical Method:
%     Second-Order Taylor Series Method
%     with Variable Mesh
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
  disp('Variable Mesh Numerical Solution for IVP Using Taylor Series Method');
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
        EU(L)= 1.0/(1+T(L)*T(L));  
    end
    % Initial conditions
 %   U(1)=EU(1);
    U(1)=1.0;
     tic  % CPU time starts
    % Solution at advanced time
       for L=1:N
           H0=H(L+1);
           T0=T(L); U0=U(L);
           DU1= -2*T0*U0*U0;
           DU2= -2*(U0*U0+2*T0*U0*DU1);
           DU3= -4*(2*U0*DU1+T0*(DU1*DU1+U0*DU2));
           DU4= -4*(3*DU1*DU1+3*U0*DU2+T0*(3*DU1*DU2+U0*DU3));
      %--- Second order approximation     
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2;
      %-----Third order approximation
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*DU3/6;
      %-----fourth order approximation
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*(DU3/6)+H0*H0*H0*H0*(DU4/24);
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
  figure (1);
  plot(T,U,'-g',T,EU,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title ('Exact vs Numerical Solution')
