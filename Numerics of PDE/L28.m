%-------------------------------- 
% Numerical Solution of 1D Wave Equation   
% Three-Level Explicit & Implicit Scheme 
% Utt=Uxx, % U(x,t)=sin(pi*x)*cos(pi*t) 
% Name:- Kshitij Bahukhandi; Roll No:-06  
% APM (III-Semester)-2024 
% SOUTH  ASIAN  UNIVERSITY 
% MATLAB File:- 28 
%----------------------------- 
  clear all;  % clear all variables in memory 
  P=0.8; 
  P0=(1-P*P)/12; 
  LL=8;   % LL=16, LL=32, LL=64 
  JJ=10;  % JJ=20, JJ=40, JJ=80,  %Solution at t=1  
  H=1.0/LL; 
  K=P*H; 
  disp('Numerical Solution of 1D Wave Equation'); 
  disp('    '); 
  fprintf('LL=%4.2f, K=%6.4f, P=%6.2f \n',LL,K,P); 
  disp('  '); 
  % Dimensions 
  X=zeros(1,LL+1);T=zeros(1,JJ+1);NU=zeros(LL+1,JJ+1); 
  
SUB=zeros(1,LL+1);DIAG=zeros(1,LL+1);SUP=zeros(1,LL+1);RH=zeros(1,LL+1); 
  U=zeros(LL+1,JJ+1);EU=zeros(LL+1,JJ+1);  
  % Formation of Nodal Points 
    for L=1:LL+1 
        for J=1:JJ+1 
                X(L)=(L-1)*H; 
                T(J)=(J-1)*K; 
        end     
    end  
    % Exact Solution & RHS function 
    for L=1:LL+1 
        for J=1:JJ+1 
                 EU(L,J)=cos(pi*T(J))*sin(pi*X(L)); 
        end 
    end 
    % Initial conditions 
    for L=1:LL+1 
        U(L,1)=EU(L,1); 
        U(L,2)=EU(L,2); 
        NU(L,1)=EU(L,1); 
        NU(L,2)=EU(L,2); 
    end
    %  Solution at advanced time level 
     tic  % CPU time starts 
    for J=2:JJ 
     % Boundary conditions 
        U(1,3)=EU(1,J+1); 
        U(LL+1,3)=EU(LL+1,J+1);  
       % Solution at Internal grid points at each time level 
       for L=2:1:LL 
           R0=U(L,2); 
           R1=U(L+1,2); 
           R2=U(L-1,2); 
           R3=R1+R2-2*R0; 
           S0=U(L,1); 
           S1=U(L+1,1); 
           S2=U(L-1,1); 
           S3=S1+S2-2*S0; 
           SUB(L)=P0; 
           DIAG(L)=1-2*P0; 
           SUP(L)=P0; 
           RH(L)=P*P*R3+2*R0+2*P0*R3-S0-P0*S3;  
       end 
           RH(2)=RH(2)-SUB(2)*U(1,3); 
           RH(LL)=RH(LL)-SUP(LL)*U(LL+1,3); 
           SUB(2)=0.0; 
           SUP(LL)=0.0; 
       % Tri-diagonal solver 
           SUP(2)=SUP(2)/DIAG(2); 
           RH(2)=RH(2)/DIAG(2); 
       for L=3:LL-1 
           SUP(L)=SUP(L)/(DIAG(L)-SUB(L)*SUP(L-1)); 
       end 
       for L=3:LL 
           RH(L)=(RH(L)-SUB(L)*RH(L-1))/(DIAG(L)-SUB(L)*SUP(L-1)); 
       end 
           U(LL,3)=RH(LL); 
       for L=2:LL-1 
           U(LL-L+1,3)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+2,3); 
       end 
       for L=1:LL+1 
           NU(L,J+1)=U(L,3); 
       end 
  
  % Calculate Errors 
           MAU=0.0; 
           RMU=0.0; 
        for L=2:LL 
            ERU=abs(NU(L,J+1)-EU(L,J+1)); 
            RMU=RMU+ERU^2; 
            if MAU<=ERU, MAU=ERU; end 
        end 
             RMU=sqrt(RMU/(LL-1)); 
     % Print errors 
       fprintf('T=%6.4f,MAU=%10.4e, RMU=%10.4e\n',J*K,MAU,RMU); 
      % Set New Initial conditions 
        for L=1:LL+1 
            U(L,1)=U(L,2); 
            U(L,2)=U(L,3); 
        end 
    end 
  disp('  '); 
  toc  % CPU time ends 
  % Plot the results 
  figure (1); 
colormap spring; 
surf(T,X,NU) 
xlabel('T values'); 
ylabel('X values'); 
zlabel('NU values'); 
title('Numerical Solution') 
  figure (2); 
colormap autumn; 
surf(T,X,EU) 
xlabel('T values'); 
ylabel('X Values'); 
zlabel('EU values'); 
title('Exact Solution') 
  figure(3); 
  plot(X,NU,'--g');  
  title('Numerical Solution') 
  figure (4); 
  plot(X,EU,'--r'); 
  title ('Exact Solution') 
  % Plot at t=1 
  for L=1:LL+1 
      US(L)=NU(L,JJ+1); 
      EUS(L)=exp(-1)*sin(pi*X(L)); 
  end
  figure (5); 
plot(X,US,'-g',X,EUS, '--r'); 
xlabel('X values'); 
ylabel('EUS values'); 
title('Numerical Vs Exact Solution at t=1') 