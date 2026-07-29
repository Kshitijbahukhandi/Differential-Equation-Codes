%--------------------------------
% Numerical Solution of 2D Heat Conduction Equation with Forcing Function
% ADI Method
% Uxx+Uyy=Ut+f(x,y,t), % U(x,y,t)=sinh(t)*sin(pi*x)*sin(pi*y)
% Name:- KSHITIJ BAHUKHANDI
% Enrollment No :- SAU/AM(M)/2023/06
% (III-Semester)-2024
% SOUTH  ASIAN  UNIVERSITY
% MATLAB File:- 15
%-----------------------------
  clear all;  % clear all variables in memory
  LAM=1.6;
  LL=8; MM=LL;  % LL=16, LL=32, Ll=64
  JJ=40;  % JJ=160, JJ=640, JJ=2560,  %Solution at t=1 
  H=1.0/LL;
  K=LAM*H*H;
  disp('Numerical Solution of 2D Heat Conduction Equation');
  disp('    ');
  fprintf('LL=%4.2f, K=%6.4f, LAM=%6.2f \n',LL,K,LAM);
  disp('  ');
  % Dimensions
  X=zeros(1,LL+1);Y=zeros(1,MM+1); T=zeros(1,JJ+1);TT=zeros(1,JJ+1);
  SUB=zeros(1,LL+1);DIAG=zeros(1,LL+1);SUP=zeros(1,LL+1);RH=zeros(1,LL+1);
  U=zeros(LL+1,MM+1,JJ+1);EU=zeros(LL+1,MM+1,JJ+1);US=zeros(LL+1,MM+1);
  % Formation of Nodal Points
    for L=1:LL+1
        for M=1:MM+1
            for J=1:JJ+1
                X(L)=(L-1)*H;
                Y(M)=(M-1)*H;
                T(J)=(J-1)*K;
                TT(J)=T(J)+0.5*K;
            end
        end      
    end 
    % Exact Solution & RHS function
    for L=1:LL+1
        for M=1:MM+1
            for J=1:JJ+1
                EU(L,M,J)=sinh(T(J))*sin(pi*X(L))*sin(pi*Y(M));
               EUT(L,M,J)=cosh(T(J))*sin(pi*X(L))*sin(pi*Y(M)); 
               EEU(L,M,J)=sinh(TT(J))*sin(pi*X(L))*sin(pi*Y(M)); 
              EEUT(L,M,J)=cosh(TT(J))*sin(pi*X(L))*sin(pi*Y(M));
                 F(L,M,J)= -2*pi*pi*EEU(L,M,J)-EEUT(L,M,J);
            end
        end
    end
    % Initial conditions
    for L=1:LL+1
        for M=1:MM+1
            U(L,M,1)=EU(L,M,1);
        end
    end
    % Boundary conditions
    for J=1:JJ
        for L=1:LL+1
            for M=1:MM+1
                U(1,M,J+1)=EU(1,M,J+1);
                U(LL+1,M,J+1)=EU(LL+1,M,J+1); 
                U(L,1,J+1)=EU(L,1,J+1);
                U(L,MM+1,J+1)=EU(L,MM+1,J+1);
            end
        end
    end
    %  Solution at advanced time level
     tic  % CPU time starts
       % LEES ADI Method
       for J=1:JJ
       % Intermediate Boundary Conditions
       for L=2:LL
           for M=1:MM:MM+1   %IB=1,2
     US(L,M)=-0.5*LAM*U(L-1,M,J+1)+(1+LAM)*U(L,M,J+1)-0.5*LAM*U(L+1,M,J+1);
           end
       end 
       % Solution at Intermediate Step
       for L=2:LL
           for M=2:MM
               F0=F(L,M,J);
               R0=U(L,M,J);
               R1=U(L+1,M,J);
               R2=U(L-1,M,J);
               R3=U(L,M+1,J);
               R4=U(L,M-1,J);
               R5=U(L+1,M+1,J);
               R6=U(L+1,M-1,J);
               R7=U(L-1,M+1,J);
               R8=U(L-1,M-1,J);
               R10=R1+R2+R3+R4;
               R11=R1+R2-2*R0;
               R12=R3+R4-2*R0;
               R18=R5+R6+R7+R8-2*R10+4*R0;
               SUB(M)=-0.5*LAM;
               DIAG(M)=1+LAM;
               SUP(M)=-0.5*LAM;
               RH(M)=R0+0.5*LAM*(R11+R12)+0.25*LAM*LAM*R18-K*F0;
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
               SUB(L)=-0.5*LAM;
               DIAG(L)=1+LAM;
               SUP(L)=-0.5*LAM;
               RH(L)=US(L,M); 
           end
               RH(2)=RH(2)-SUB(2)*U(1,M,J+1);
               RH(LL)=RH(LL)-SUP(LL)*U(LL+1,M,J+1);
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
               U(LL,M,J+1)=RH(LL);
           for L=2:LL-1
               U(LL-L+1,M,J+1)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+2,M,J+1);
           end
       end
  % Calculate Errors
               MAU=0.0;
               RMU=0.0;
        for L=2:LL
            for M=2:MM
                ERU=abs(U(L,M,J+1)-EU(L,M,J+1));
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
              U(L,M,J)=U(L,M,J+1);
          end
      end
    end
  disp('  ');
   toc  % CPU time ends
     % Plot the results at t=1
  for L=1:LL+1
      for M=1:MM+1
          UU(L,M)=U(L,M,JJ+1);
         EUU(L,M)=sinh(1)*sin(pi*X(L))*sin(pi*Y(M));
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