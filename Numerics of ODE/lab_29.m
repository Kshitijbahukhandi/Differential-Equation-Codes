%--------------------------------
% Name:  KSHITIJ BAHUKHANDI
% Roll No : SAU/AM(M)/2023/06
% M.Sc AM (II-Semester)-2024
% SOUTH  ASIAN  UNIVERSITY
% MATLAB File:- 29
% Numerical Solution ofsystem of IVPs - Taylor Series Method
% Ut= -3u+2v ;  u=(exp(-t)-exp(-6t))/5
% Vt= 3u-4v  ;  v=(2*exp(-t)+3*exp(-6t))/10
%-----------------------------


  clear all;  % clear all variables in memory
  N=10; 
  for I=1:8
  H=1.0/N;
  disp('Numerical Solution for system of IVPs Using Taylor Series Method');
  disp('    ');
  fprintf('N=%4.2f, H=%6.4f\n',N,H);
  disp('  ');
  % Formation of Nodal Points
    for L=1:N+1
          T(L)=(L-1)*H;   
    end 
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
           T0=T(L); U0=U(L); V0=V(L);
           DU1= -3*U0+2*V0;
           DV1=  3*U0-4*V0;
           DU2= -3*DU1+2*DV1;
           DV2= 3*DU1-4*DV1;
           DU3= -3*DU2+2*DV2;
           DV3= 3*DU2-4*DV2;
           DU4= -3*DU3+2*DV3;
           DV4= 3*DU3-4*DV3;
      %--- Second order approximation     
           U(L+1)=U(L)+H*DU1+0.5*H*H*DU2;
           V(L+1)=V(L)+H*DV1+0.5*H*H*DV2;
      %-----Third order approximation
           U(L+1)=U(L)+H*DU1+0.5*H*H*DU2+H*H*H*DU3/6;
           V(L+1)=V(L)+H*DV1+0.5*H*H*DV2+H*H*H*DV3/6;
      %-----fourth order approximation
           U(L+1)=U(L)+H*DU1+0.5*H*H*DU2+H*H*H*(DU3/6)+H*H*H*H*(DU4/24);
           V(L+1)=V(L)+H*DV1+0.5*H*H*DV2+H*H*H*(DV3/6)+H*H*H*H*(DV4/24);
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
