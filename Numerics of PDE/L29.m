%-------------------------------- 
% Numerical Solution of 1D Nonlinear Hyperbolic Eqn with a forcing function  
% Utt=Uxx-alpha*Ut*Ut*Ut+f(x,t), % U(x,t)=sinh(t)*sin(pi*x) 
% Name:- Kshitij Bahukhandi 
% Roll No:- SAU/AM(M)/2023/06 
% APM (III-Semester)- 2024 
% SOUTH  ASIAN  UNIVERSITY 
% MATLAB File:- 29 
%----------------------------- 
  clear all;  % clear all variables in memory 
  ALPHA=10.0; 
  P=0.8; 
  LL=8;   % LL=16, LL=32, Ll=64 
  JJ=10;  % JJ=20, JJ=40, JJ=80,  %Solution at t=1.0  
  H=1.0/LL;   
  K=P*H; 
  disp('Numerical Solution of 1D Nonlinear Hyperbolic Equation'); 
  disp('    '); 
  fprintf('LL=%4.2f, K=%6.4f, P=%6.2f \n',LL,K,P); 
  disp('  '); 
  % Dimensions 
  
X=zeros(1,LL+1);T=zeros(1,JJ+1);NU=zeros(LL+1,JJ+1);F=zeros(LL+1,JJ+1); 
  U=zeros(LL+1,JJ+1);EU=zeros(LL+1,JJ+1);OLDU=zeros(LL+1,JJ+1); 
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
            EU(L,J)=sinh(T(J))*sin(pi*X(L)); 
           EUT(L,J)=cosh(T(J))*sin(pi*X(L)); 
            F(L,J)=(1+pi*pi)*EU(L,J)+ALPHA*EUT(L,J)*EUT(L,J); 
        end 
    end 
    % Initial conditions 
    for L=1:LL+1 
        U(L,1)=EU(L,1); 
        U(L,2)=EU(L,2); 
        NU(L,1)=EU(L,1); 
        NU(L,2)=EU(L,2); 
    end 
   % Boundary conditions 
     for J=3:JJ+1 
         U(1,J)=EU(1,J); 
         U(LL+1,J)=EU(LL+1,J);  
     end 
    %  Solution at advanced time level 
     tic  % CPU time starts 
    for J=2:JJ 
    % Initial & First Approx 
    for L=2:LL 
        OLDU(L,3)=EU(L,J+1); 
        U(L,3)=OLDU(L,3); 
    end 
    for KK=1:6000   %Iteration Starts 
    % Solution at Internal grid points at each time level 
       for L=2:LL 
           F0=F(L,J); 
           U0=U(L,3); 
           U1=U(L+1,3); 
           U2=U(L-1,3); 
           R0=U(L,2); 
           R1=U(L+1,2); 
           R2=U(L-1,2); 
           R3=R1+R2-2*R0; 
           S0=U(L,1); 
           S1=U(L+1,1); 
           S2=U(L-1,1); 
           AA=0.25*ALPHA*U0*U0+(1-0.5*ALPHA*S0)*U0+S0-2*R0; 
           AA=AA+0.25*ALPHA*S0*S0-P*P*R3-K*K*F0; 
           DA=0.5*ALPHA*U0+1-0.5*ALPHA*S0; 
           U(L,3)=U(L,3)-AA/DA;  
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
      % fprintf('T=%6.4f,MAU=%10.4e, RMU=%10.4e\n',(J-1)*K,MAU,RMU);
      % Error Tolerance 
            BIG=0.0; 
        for L=2:LL 
            FU=abs(U(L,3)-OLDU(L,3)); 
            if FU >= BIG, BIG=FU; end 
        end 
            if BIG <= 10^(-12), break; end 
        for L=2:LL 
            OLDU(L,3)=U(L,3); 
        end  
    end 
      % Print errors 
       fprintf('T=%6.4f,KK=%4.0f,MAU=%10.4e, RMU=%10.4e\n',J*K,KK,MAU,RMU); 
      % Set New Initial Conditions 
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
      EUS(L)=sinh(1)*sin(pi*X(L)); 
  end 
  figure (5); 
  plot(X,US,'-g',X,EUS, '--r'); 
  xlabel('X values'); 
  ylabel('EUS values'); 
  title('Numerical Vs Exact Solution at t=1') 