%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Taylor Series Method on a Variable Mesh
%
% System of Differential Equations:
%     du/dt = -3u + 2v
%     dv/dt =  3u - 4v
%
% Exact Solutions:
%     u(t) = (exp(-t) - exp(-6t))/5
%     v(t) = (2*exp(-t) + 3*exp(-6t))/10
%
% Numerical Method:
%     Taylor Series Method with Variable Mesh
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
  H=1.0/N;
  disp('Variable Mesh for system of IVPs Using Taylor Series Method');
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
        EU(L)= (exp(-T(L))-exp(-6*T(L)))/5; 
        EV(L)= (2*exp(-T(L))+3*exp(-6*T(L)))/10; 
    end
    % Initial conditions
 %   U(1)=EU(1);
 %   V(1)=EV(1);
    U(1)=0.0;
    V(1)=0.5;
     tic  % CPU time starts
    % Solution at advanced time
       for L=1:N 
           T0=T(L); U0=U(L); V0=V(L); H0=H(L+1);
           DU1= -3*U0+2*V0;
           DV1=  3*U0-4*V0;
           DU2= -3*DU1+2*DV1;
           DV2= 3*DU1-4*DV1;
           DU3= -3*DU2+2*DV2;
           DV3= 3*DU2-4*DV2;
           DU4= -3*DU3+2*DV3;
           DV4= 3*DU3-4*DV3;
      %--- Second order approximation     
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2;
           V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2;
      %-----Third order approximation
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*DU3/6;
           V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2+H0*H0*H0*DV3/6;
      %-----fourth order approximation
           U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*(DU3/6)+H0*H0*H0*H0*(DU4/24);
           V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2+H0*H0*H0*(DV3/6)+H0*H0*H0*H0*(DV4/24);
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

        
        
