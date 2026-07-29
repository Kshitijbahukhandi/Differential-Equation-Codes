        %--------------------------------
% Name:  KSHITIJ BAHUKHANDI
% Roll No : SAU/AM(M)/2023/06
% M.Sc AM (II-Semester)-2024
% SOUTH  ASIAN  UNIVERSITY
% MATLAB File:- 28
% Numerical Solution of IVP on a Variable Mesh - 4th order explicit RK-method 
% Ut= -2*t*U*U
% U(t)= 1/(1+t*t)
%-----------------------------


  clear all;  % clear all variables in memory
  N=10; 
  ETA=1.2;
  disp('Variable Mesh Numerical Solution for IVP Using 4th order explicit RK Method');
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
      %---Classical Method
           K1= -2*H0*T0*U0*U0;
           K2= -2*H0*(T0+H0/2)*(U0+K1/2)*(U0+K1/2);
           K3= -2*H0*(T0+H0/2)*(U0+K2/2)*(U0+K2/2);
           K4= -2*H0*(T0+H0)*(U0+K3)*(U0+K3);
           U(L+1)=U(L)+(K1+2*K2+2*K3+K4)/6;    
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
