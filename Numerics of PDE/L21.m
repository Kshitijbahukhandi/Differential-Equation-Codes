%--------------------------------
%Name= Kshitij Bahukhandi
%Roll NO.= SAU/AM(M)/2022/06
%LAB:- 21
%SOUTH ASIAN UNIVERSITY
% Numerical Solution of 2D Wave Equation - ADI Method of O(k2+h2)
% Unconditionally Stable
% Utt=Uxx+Uyy, % U(x,y,t)=cos(sqrt(2)*pi*t)*sin(pi*x)*sin(pi*y)
%-----------------------------
  clear all;  % clear all variables in memory
  P=1.6;
  GAMA=0.5;
  P0=-GAMA*P*P;
  LL=32;   % LL=16, LL=32, Ll=64
  MM=LL;
  JJ=20;  % JJ=20, JJ=40, JJ=80,  %Solution at t=1 
  H=1.0/LL;
  K=P*H;
  disp('Numerical Solution of 2D Wave Equation - ADI Method');
  disp('    ');
  fprintf('LL=%4.2f, K=%6.4f, P=%6.2f \n',LL,K,P);
  disp('  ');
  % Dimensions
  X=zeros(1,LL+1);Y=zeros(1,MM+1); T=zeros(1,JJ+1); 
  SUB=zeros(1,LL+1);DIAG=zeros(1,LL+1);SUP=zeros(1,LL+1);RH=zeros(1,LL+1);
  U=zeros(LL+1,MM+1,JJ+1);EU=zeros(LL+1,MM+1,JJ+1);US=zeros(LL+1,MM+1);
  NU=zeros(LL+1,MM+1,JJ+1);
  % Formation of Nodal Points
    for L=1:LL+1
        for M=1:MM+1
            for J=1:JJ+1
                X(L)=(L-1)*H;
                Y(M)=(M-1)*H;
                T(J)=(J-1)*K; 
            end
        end      
    end 
    % Exact Solution & RHS function
    for L=1:LL+1
        for M=1:MM+1
            for J=1:JJ+1
                EU(L,M,J)=cos(sqrt(2)*pi*T(J))*sin(pi*X(L))*sin(pi*Y(M));
            end
        end
    end
    % Initial conditions
    for L=1:LL+1
        for M=1:MM+1
            U(L,M,1)=EU(L,M,1);
            U(L,M,2)=EU(L,M,2);
           NU(L,M,1)=EU(L,M,1);
           NU(L,M,2)=EU(L,M,2);
        end
    end
    %  Solution at advanced time level
     tic  % CPU time starts
    for J=2:JJ
       % LEES ADI Method
         % Boundary conditions
         for L=1:LL+1
             for M=1:MM+1
                 U(L,1,3)=EU(L,1,J+1);
                 U(L,MM+1,3)=EU(L,MM+1,J+1);
                 U(1,M,3)=EU(1,M,J+1);
                 U(LL+1,M,3)=EU(LL+1,M,J+1);
             end
         end
         % Intermediate Boundary Conditions
       for L=2:LL
           for M=1:MM:MM+1   %IB=1,2
               U0=U(L,M,3);
               U1=U(L+1,M,3);
               U2=U(L-1,M,3);
               R0=U(L,M,2);
               R1=U(L+1,M,2);
               R2=U(L-1,M,2);
               S0=U(L,M,1);
               S1=U(L+1,M,1);
               S2=U(L-1,M,1);
               T0=U0-2*R0+S0;
               U9=U1+U2-2*U0;
               R9=R1+R2-2*R0;
               S9=S1+S2-2*S0;
               US(L,M)=T0+P0*(U9-2*R9+S9);
           end
       end 
       % Solution at Intermediate Step
       for L=2:LL
           for M=2:MM
               %F0=F(L,M,J);
               R0=U(L,M,2);
               R1=U(L+1,M,2);
               R2=U(L-1,M,2);
               R3=U(L,M+1,2);
               R4=U(L,M-1,2);
               R9=R1+R2-2*R0;
               R10=R3+R4-2*R0;
               SUB(M)=P0;
               DIAG(M)=1-2*P0;
               SUP(M)=P0;
               RH(M)=P*P*(R9+R10);
           end     
               RH(2)=RH(2)-SUB(2)*US(L,1);
               RH(MM)=RH(MM)-SUP(MM)*US(L,MM+1);
               SUB(2)=0.0;
               SUP(MM)=0.0;
       % Tri-diagonal solver
               SUP(2)=SUP(2)/DIAG(2);
               RH(2)=RH(2)/DIAG(2);
           for M=3:MM-1
               SUP(M)=SUP(M)/(DIAG(M)-SUB(M)*SUP(M-1));
           end
           for M=3:MM
               RH(M)=(RH(M)-SUB(M)*RH(M-1))/(DIAG(M)-SUB(M)*SUP(M-1));
           end
               US(L,MM)=RH(MM);
           for M=2:MM-1
               US(L,MM-M+1)=RH(MM-M+1)-SUP(MM-M+1)*US(L,MM-M+2);
           end
       end
       % Solution at Main Time Step
       for M=2:MM
           for L=2:LL
               R0=U(L,M,2);
               R1=U(L+1,M,2);
               R2=U(L-1,M,2);
               S0=U(L,M,1);
               S1=U(L+1,M,1);
               S2=U(L-1,M,1);
               R9=R1+R2-2*R0;
               S9=S1+S2-2*S0;
               SUB(L)=P0;
               DIAG(L)=1-2*P0;
               SUP(L)=P0;
               RH(L)=US(L,M)+2*R0-S0+P0*(2*R9-S9); 
           end
               RH(2)=RH(2)-SUB(2)*U(1,M,3);
               RH(LL)=RH(LL)-SUP(LL)*U(LL+1,M,3);
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
               U(LL,M,3)=RH(LL);
           for L=2:LL-1
               U(LL-L+1,M,3)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+2,M,3);
           end
       end
       %  Numerical solution for Graph
       for L=1:LL+1
           for M=1:MM+1
               NU(L,M,J+1)=U(L,M,3);
           end
       end
  % Calculate Errors
               MAU=0.0;
               RMU=0.0;
        for L=2:LL
            for M=2:MM
                ERU=abs(NU(L,M,J+1)-EU(L,M,J+1));
                RMU=RMU+ERU^2;
                if MAU<=ERU, MAU=ERU; end
            end
        end
                RMU=sqrt(RMU/((LL-1)*(MM-1)));
     % Print errors
       fprintf('T=%6.4f,MAU=%10.4e, RMU=%10.4e\n',J*K,MAU,RMU);
      % Set New Initial Conditions
      for L=1:LL+1
          for M=1:MM+1
              U(L,M,1)=U(L,M,2);
              U(L,M,2)=U(L,M,3);
          end
      end
    end
  disp('  ');
   toc  % CPU time ends
     % Plot the results at t=1
  for L=1:LL+1
      for M=1:MM+1
          UU(L,M)=NU(L,M,JJ+1);
          EUU(L,M)=cos(sqrt(2)*pi)*sin(pi*X(L))*sin(pi*Y(M));
      end
  end
  figure (1);
  colormap summer;
  surf(X,Y,UU);
  xlabel('X values');
  ylabel('Y values');
  zlabel('UU values');
  title('Numerical Solution at t=1')
  figure (2);
  colormap winter;
  surf(X,Y,EUU);
  xlabel('X values');
  ylabel('Y values');
  zlabel('EUU values');
  title('Exact Solution at t=1')


         
       
        