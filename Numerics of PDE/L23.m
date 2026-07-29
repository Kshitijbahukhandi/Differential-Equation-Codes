%--------------------------------
%Name= Kshitij Bahukhandi
%Roll NO.= SAU/AM(M)/2022/06
%LAB:- 23
%SOUTH ASIAN UNIVERSITY
% Numerical Solution of 2D Nonlinear Hyperbolic Eqn with a forcing function 
% Utt=Uxx+Uyy-alpha*Ut*Ut+f(x,y,t), % U(x,y,t)=sinh(t)*sin(pi*x)*sin(pi*y)
%-----------------------------
  clear all;  % clear all variables in memory
  ALPHA=10.0;
  P=0.8;
  LL=16;   % LL=16, LL=32, Ll=64
  MM=LL;
  JJ=20;  % JJ=20, JJ=40, JJ=80,  %Solution at t=1.0 
  H=1.0/LL;
  K=P*H;
  disp('Numerical Solution of 2D Nonlinear Hyperbolic Equation');
  disp('    ');
  fprintf('LL=%4.2f, K=%6.4f, P=%6.2f \n',LL,K,P);
  disp('  ');
  % Dimensions
  X=zeros(1,LL+1); Y=zeros(1,MM+1); T=zeros(1,JJ+1);
  U=zeros(LL+1,MM+1,JJ+1);OLDU=zeros(LL+1,MM+1,JJ+1);
  NU=zeros(LL+1,MM+1,JJ+1);EU=zeros(LL+1,MM+1,JJ+1);
  F=zeros(LL+1,MM+1,JJ+1);
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
            EU(L,M,J)=sinh(T(J))*sin(pi*X(L))*sin(pi*Y(M));
           EUT(L,M,J)=cosh(T(J))*sin(pi*X(L))*sin(pi*Y(M));
            F(L,M,J)=(1+2*pi*pi)*EU(L,M,J)+ALPHA*EUT(L,M,J)*EUT(L,M,J);
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
   % Boundary conditions
     for J=3:JJ+1
         for L=1:LL+1
             for M=1:MM+1
                 U(L,1,J)=EU(L,1,J);
                 U(L,MM+1,J)=EU(L,MM+1,J);
                 U(1,M,J)=EU(1,M,J);
                 U(LL+1,M,J)=EU(LL+1,M,J);
             end
         end
     end
    %  Solution at advanced time level
     tic  % CPU time starts
    for J=2:JJ
    % Initial & First Approx
    for L=2:LL
        for M=2:MM
            OLDU(L,M,3)=EU(L,M,J+1);
            U(L,M,3)=OLDU(L,M,3);
        end
    end
    for KK=1:600   %Iteration Starts
    % Solution at Internal grid points at each time level
       for L=2:LL
           for M=2:MM
               F0=F(L,M,J);
               U0=U(L,M,3);
               U1=U(L+1,M,3);
               U2=U(L-1,M,3);
               U3=U(L,M+1,3);
               U4=U(L,M-1,3);
               R0=U(L,M,2);
               R1=U(L+1,M,2);
               R2=U(L-1,M,2);
               R3=U(L,M+1,2);
               R4=U(L,M-1,2);
               R5=R1+R2-2*R0;
               R6=R3+R4-2*R0;
               S0=U(L,M,1);
               S1=U(L+1,M,1);
               S2=U(L-1,M,1);
               S3=U(L,M+1,1);
               S4=U(L,M-1,1);
               AA=0.25*ALPHA*U0*U0+(1-0.5*ALPHA*S0)*U0+S0-2*R0;
               AA=AA+0.25*ALPHA*S0*S0-P*P*(R5+R6)-K*K*F0;
               DA=0.5*ALPHA*U0+1-0.5*ALPHA*S0;
               U(L,M,3)=U(L,M,3)-AA/DA; 
           end
       end
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
      % fprintf('T=%6.4f,MAU=%10.4e, RMU=%10.4e\n',J*K,MAU,RMU);
      % Error Tolerance
            BIG=0.0;
        for L=2:LL
            for M=2:MM
                FU=abs(U(L,M,3)-OLDU(L,M,3));
            if FU >= BIG, BIG=FU; end
            end
        end
            if BIG <= 10^(-12), break; end
        for L=2:LL
            for M=2:MM
                OLDU(L,M,3)=U(L,M,3);
            end
        end 
    end
      % Print errors
       fprintf('T=%6.4f,KK=%4.0f,MAU=%10.4e, RMU=%10.4e\n',J*K,KK,MAU,RMU);
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
   
        