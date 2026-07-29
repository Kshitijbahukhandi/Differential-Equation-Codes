%--------------------------------
%Name= Kshitij Bahukhandi
%Roll NO.= SAU/AM(M)/2022/06
%LAB:- 25
%SOUTH ASIAN UNIVERSITY
% Numerical Solution of 3D Wave Equation - ADI Method of O(k4+h4)
% Utt=Uxx+Uyy+Uzz %U=cos(sqrt(3)*pi*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
%-----------------------------
  clear all;  % clear all variables in memory
  P=0.4;
  P0=(1-P*P)/12;
  LL=16;          % LL=16, LL=32, Ll=64
  MM=LL;
  NN=LL;
  JJ=40;  % JJ=20, JJ=40, JJ=80,  %Solution at t=1 
  H=1.0/LL;
  K=P*H;
  disp('Numerical Solution of 3D Wave Equation - ADI Method');
  disp('    ');
  fprintf('LL=%4.2f, K=%6.4f, P=%6.2f \n',LL,K,P);
  disp('  ');
  % Dimensions
  X=zeros(1,LL+1);Y=zeros(1,MM+1);Z=zeros(1,NN+1); T=zeros(1,JJ+1);
  SUB=zeros(1,LL+1);DIAG=zeros(1,LL+1);SUP=zeros(1,LL+1);RH=zeros(1,LL+1);
  U=zeros(LL+1,MM+1,NN+1,JJ+1);EU=zeros(LL+1,MM+1,NN+1,JJ+1);
  US=zeros(LL+1,MM+1,NN+1);USS=zeros(LL+1,MM+1,NN+1);
  F=zeros(LL+1,MM+1,NN+1,JJ+1);NU=zeros(LL+1,MM+1,NN+1,JJ+1);
  % Formation of Nodal Points
    for L=1:LL+1
        for M=1:MM+1
            for N=1:NN+1
                for J=1:JJ+1
                    X(L)=(L-1)*H;
                    Y(M)=(M-1)*H;
                    Z(N)=(N-1)*H;
                    T(J)=(J-1)*K;
                end
            end
        end      
    end 
    % Exact Solution & RHS function
    for L=1:LL+1
        for M=1:MM+1
            for N=1:NN+1
                for J=1:JJ+1
   EU(L,M,N,J)=cos(sqrt(3)*pi*T(J))*sin(pi*X(L))*sin(pi*Y(M))*sin(pi*Z(N));
                end
            end
        end
    end
    % Initial conditions
    for L=1:LL+1
        for M=1:MM+1
            for N=1:NN+1
                U(L,M,N,1)=EU(L,M,N,1);
                U(L,M,N,2)=EU(L,M,N,2);
               NU(L,M,N,1)=EU(L,M,N,1);
               NU(L,M,N,2)=EU(L,M,N,2);
            end
        end
    end
    %  Solution at advanced time level
     tic  % CPU time starts
    for J=2:JJ
       % Lee's ADI Method
          % Boundary conditions
        for L=1:LL+1
            for M=1:MM+1
                for N=1:NN+1
                    U(1,M,N,3)=EU(1,M,N,J+1);
                    U(LL+1,M,N,3)=EU(LL+1,M,N,J+1); 
                    U(L,1,N,3)=EU(L,1,N,J+1);
                    U(L,MM+1,N,3)=EU(L,MM+1,N,J+1);
                    U(L,M,1,3)=EU(L,M,1,J+1);
                    U(L,M,NN+1,3)=EU(L,M,NN+1,J+1);
                end
            end
        end
       % Second Intermediate Boundary Conditions for US
       for L=2:LL 
           for M=1:MM:MM+1
               for N=2:NN
                   U0=U(L,M,N,3);
                   U1=U(L+1,M,N,3);
                   U2=U(L-1,M,N,3);
                   R0=U(L,M,N,2);
                   R1=U(L+1,M,N,2);
                   R2=U(L-1,M,N,2);
                   S0=U(L,M,N,1);
                   S1=U(L+1,M,N,1);
                   S2=U(L-1,M,N,1);
                   T20=U0-2*R0+S0;
                   U20=U1+U2-2*U0;
                   R20=R1+R2-2*R0;
                   S20=S1+S2-2*S0;
                   US(L,M,N)=T20+P0*(U20-2*R20+S20);
               end
           end
       end 
       % First Intermediate Boundary Conditions for USS
       for L=2:LL
           for M=2:MM   
           for N=1:NN:NN+1
               USS(L,M,N)=P0*US(L,M-1,N)+(1-2*P0)*US(L,M,N)+P0*US(L,M+1,N);
           end
           end
       end 
       % Solution at First Intermediate Step
       for L=2:LL
           for M=2:MM
               for N=2:NN
                   R0=U(L,M,N,2);
                   R1=U(L+1,M,N,2);
                   R2=U(L-1,M,N,2);
                   R3=U(L,M+1,N,2);
                   R4=U(L,M-1,N,2);
                   R5=U(L,M,N+1,2);
                   R6=U(L,M,N-1,2);
                   R7=U(L+1,M+1,N,2);
                   R8=U(L+1,M-1,N,2);
                   R9=U(L-1,M+1,N,2);
                   R10=U(L-1,M-1,N,2);
                   R11=U(L,M+1,N+1,2);
                   R12=U(L,M+1,N-1,2);
                   R13=U(L,M-1,N+1,2);
                   R14=U(L,M-1,N-1,2);
                   R15=U(L+1,M,N+1,2);
                   R16=U(L+1,M,N-1,2);
                   R17=U(L-1,M,N+1,2);
                   R18=U(L-1,M,N-1,2);
                   R20=R1+R2-2*R0;
                   R21=R3+R4-2*R0;
                   R22=R5+R6-2*R0;
                   R23=(R7+R8+R9+R10-2*(R1+R2+R3+R4)+4*R0)/6;
                   R24=(R11+R12+R13+R14-2*(R3+R4+R5+R6)+4*R0)/6;
                   R25=(R15+R16+R17+R18-2*(R1+R2+R5+R6)+4*R0)/6;
                   SUB(N)=P0;
                   DIAG(N)=1-2*P0;
                   SUP(N)=P0;
                   RH(N)=P*P*(R20+R21+R22+R23+R24+R25);
               end     
                   RH(2)=RH(2)-SUB(2)*USS(L,M,1);
                   RH(NN)=RH(NN)-SUP(NN)*USS(L,M,NN+1);
                   SUB(2)=0.0;
                   SUP(NN)=0.0;
       % Tri-diagonal solver
                   SUP(2)=SUP(2)/DIAG(2);
                   RH(2)=RH(2)/DIAG(2);
               for N=3:NN-1
                   SUP(N)=SUP(N)/(DIAG(N)-SUB(N)*SUP(N-1));
               end
               for N=3:NN
                   RH(N)=(RH(N)-SUB(N)*RH(N-1))/(DIAG(N)-SUB(N)*SUP(N-1));
               end
                   USS(L,M,NN)=RH(NN);
               for N=2:NN-1
                   USS(L,M,NN-N+1)=RH(NN-N+1)-SUP(NN-N+1)*USS(L,M,NN-N+2);
               end
           end
       end
       % Solution at Second Intermediate Step
       for L=2:LL
           for N=2:NN
               for M=2:MM
                   SUB(M)=P0;
                   DIAG(M)=1-2*P0;
                   SUP(M)=P0;
                   RH(M)=USS(L,M,N);
               end     
                   RH(2)=RH(2)-SUB(2)*US(L,1,N);
                   RH(MM)=RH(MM)-SUP(MM)*US(L,MM+1,N);
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
                   US(L,MM,N)=RH(MM);
               for M=2:MM-1
                   US(L,MM-M+1,N)=RH(MM-M+1)-SUP(MM-M+1)*US(L,MM-M+2,N);
               end
           end
       end
       % Solution at Main Time Step
       for M=2:MM
           for N=2:NN
               for L=2:LL
                   R0=U(L,M,N,2);
                   R1=U(L+1,M,N,2);
                   R2=U(L-1,M,N,2);
                   S0=U(L,M,N,1);
                   S1=U(L+1,M,N,1);
                   S2=U(L-1,M,N,1);
                   R20=R1+R2-2*R0;
                   S20=S1+S2-2*S0;
                   SUB(L)=P0;
                   DIAG(L)=1-2*P0;
                   SUP(L)=P0;
                   RH(L)=US(L,M,N)+2*R0-S0+P0*(2*R20-S20);
               end
                   RH(2)=RH(2)-SUB(2)*U(1,M,N,3);
                   RH(LL)=RH(LL)-SUP(LL)*U(LL+1,M,N,3);
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
                  U(LL,M,N,3)=RH(LL);
              for L=2:LL-1
                U(LL-L+1,M,N,3)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+2,M,N,3);
              end
           end
       end
       %  Num Solution
       for L=1:LL+1
           for M=1:MM+1
               for N=1:NN+1
                   NU(L,M,N,J+1)=U(L,M,N,3);
               end
           end
       end
  % Calculate Errors
               MAU=0.0;
               RMU=0.0;
        for L=2:LL
            for M=2:MM
                for N=2:NN
                    ERU=abs(NU(L,M,N,J+1)-EU(L,M,N,J+1));
                    RMU=RMU+ERU^2;
                    if MAU<=ERU, MAU=ERU; end
                end
            end
        end
                RMU=sqrt(RMU/((LL-1)*(MM-1)*(NN-1)));
     % Print errors
       fprintf('T=%6.4f,MAU=%10.4e, RMU=%10.4e\n',J*K,MAU,RMU);
      % Set New Initial Conditions
      for L=1:LL+1
          for M=1:MM+1
              for N=1:NN+1
                  U(L,M,N,1)=U(L,M,N,2);
                  U(L,M,N,2)=U(L,M,N,3);
              end
          end
      end
    end
  disp('  ');
   toc  % CPU time ends
    % Plot the results at t=1 and z=0.5
    NNN=0.5*NN+1;
  for L=1:LL+1
      for M=1:MM+1
        UU(L,M)=NU(L,M,NNN,JJ+1);
       EUU(L,M)=cos(sqrt(3)*pi)*sin(pi*X(L))*sin(pi*Y(M))*sin(pi*Z(NNN));
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

     
         
       
        