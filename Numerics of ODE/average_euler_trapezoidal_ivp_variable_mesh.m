%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Average Euler-Trapezoidal Method on a Variable Mesh
%
% Differential Equation:
%     du/dt = -2*t*u^2
%
% Exact Solution:
%     u(t) = 1/(1+t^2)
%
% Numerical Method:
%     Average Euler-Trapezoidal Method with Variable Mesh
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
  ETA=1.2;
  N=10; % No of internal grid points
  for I=1:8
  disp('Variable Mesh Method for Non-Linear IVP -Average Euler Method');
  disp('    ');
  fprintf('N=%4.2f\n',N);
  disp('  ');
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
    U(1)=EU(1);
    NU(1)=EU(1);
     tic  % CPU time starts
    % Solution at advanced time point
       for L=1:N
           OLDU(2)=EU(L+1);
           U(2)=OLDU(2);
           for KK=1:600
           H0=H(L+1); 
           U1=U(1); U2=U(2); 
           AA=H0*T(L+1)*U2*U2+U2+H0*T(L)*U1*U1-U1;
           DA=2*H0*T(L+1)*U2+1;
           U(2)=U(2)-AA/DA;
           NU(L+1)=U(2);
         %  ERU=abs(NU(L+1)-EU(L+1));
           BIG=abs(U(2)-OLDU(2));
           if BIG <=10^(-12), break; end
           OLDU(2)=U(2);
           end
         %   fprintf('N=%4.2f, ERU=%10.4e\n',N,ERU);
          %  disp('  ');
            %   Set new Initial Conditions
            U(1)=U(2);
       end
           ERR=abs(NU(N+1)-EU(N+1));
           % Print errors
           fprintf('N=%4.2f, ERR=%10.4e\n',N,ERR);
           disp('  ');  
           N=2*N;
  end
   toc  % CPU time ends
    % Plot the results
  figure(1);
  plot(T,NU,'--g'); 
  xlabel('t-values');
  ylabel('Numerical values u');
  title('Numerical Solution')
  figure (2);
  plot(T,EU,'--r');
  xlabel('t-values');
  ylabel('Exact values u');
  title ('Exact Solution')
  figure (3);
  plot(T,NU,'-g',T,EU,'--r');
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title ('Exact vs Numerical Solution')
  figure(4);
  plot(T,NU,'b+'); hold on
  plot(T,EU,'-r'); hold off
  xlabel('t-values');
  ylabel('Numerical and Exact values u');
  title('Numerical vs Exact Solution')
 
