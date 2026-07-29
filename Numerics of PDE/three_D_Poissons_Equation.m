%-------------------------------------------------------------------------- 
% NAME:- KSHITIJ BAHUKHANDI
% South Asian University 
% ENROLLMENT NUMBER:- SAU\AM(M)\2023\06
% 3rd Sem, Numerics of PDE 
% Lab File 05 
% Numerical Solution of 3D Poissons Equation using o(h^2) method 
% Uxx+Uyy+Uzz=f(x,y,z), % U(x,y,z)=exp(x*y*z) 
  
%-------------------------------------------------------------------------- 
  clear all;  % clear all variables in memory 
  LL=7; 
  MM=LL; 
  NN=LL; 
  H=1.0/(LL+1); 
  disp('Numerical Solution of 3D Poissons Equation'); 
  disp('    '); 
  fprintf('LL=%4.2f \n',LL); 
  disp('  '); 
   % Dimensions 
  % Formation of Nodal Points 
    for L=1:LL+2 
        for M=1:MM+2 
            for N=1:NN+2 
                X(L)=(L-1)*H; 
                Y(M)=(M-1)*H; 
                Z(N)=(N-1)*H; 
            end 
        end 
    end     
    % Exact Solution & RHS function 
    for L=1:LL+2 
        for M=1:MM+2 
            for N=1:NN+2 
                EU(L,M,N)=exp(X(L)*Y(M)*Z(N)); 
                F(L,M,N)=(X(L)^2+Y(M)^2+Z(N)^2)*EU(L,M,N); 
            end 
        end 
    end 
    % Boundary conditions 
    for L=1:LL+2 
        for M=1:MM+2 
            for N=1:NN+2 
                U(L,M,1)=EU(L,M,1); 
                U(L,M,NN+2)=EU(L,M,NN+2); 
                U(L,1,N)=EU(L,1,N); 
                U(L,MM+2,N)=EU(L,MM+2,N); 
                U(1,M,N)=EU(1,M,N); 
                U(LL+2,M,N)=EU(LL+2,M,N); 
           end 
        end  
 
    end 
    % Initial and First Approximations 
    for L=2:LL+1 
        for M=2:MM+1 
            for N=2:NN+1 
                OLDU(L,M,N)=0.0; 
                U(L,M,N)=OLDU(L,M,N); 
            end 
        end 
    end 
    tic  % CPU time starts 
    % Iteration starts  
    for KK=1:6000 
    %---------------- 
    for N=2:NN+1     
         for M=2:MM+1 
            for L=2:LL+1 
                F0=F(L,M,N); 
                FF=H*H*F0; 
                U0=U(L,M,N); 
                U1=U(L+1,M,N); 
                U2=U(L-1,M,N); 
                U3=U(L,M+1,N); 
                U4=U(L,M-1,N); 
                U5=U(L,M,N+1); 
                U6=U(L,M,N-1); 
                U(L,M,N)=(U1+U2+U3+U4+U5+U6-FF)/6; 
            end 
         end 
    end 
        % Calculate Errors 
        MAU=0.0; 
        RMU=0.0; 
        for L=2:LL+1 
            for M=2:MM+1 
                for N=2:NN+1 
                    ERU=abs(U(L,M,N)-EU(L,M,N)); 
                    RMU=RMU+ERU^2; 
                    if MAU<=ERU, MAU=ERU; end 
                end 
            end 
        end 
                RMU=sqrt(RMU/(LL^3)); 
     % Print errors 
      fprintf('KK=%4.0f,  MAU=%10.4e, RMU=%10.4e\n',KK,MAU,RMU); 
      % Error Tolerance 
      BIG=0.0; 
      for L=2:LL+1 
            for M=2:MM+1 
                for N=2:NN+1 
                    FU=abs(U(L,M,N)-OLDU(L,M,N)); 
                    if FU >= BIG, BIG=FU; end 
 
                end 
            end 
      end 
            if BIG <= 10^(-12), break; end 
       for L=2:LL+1 
             for M=2:MM+1 
                  for N=2:NN+1 
                      OLDU(L,M,N)=U(L,M,N); 
                  end 
             end 
       end 
    end 
 disp('  '); 
 toc  % CPU time ends 
  % Plot the results at z=0.5 
    NNN=0.5*(NN+1); 
  for L=1:LL+2 
      for M=1:MM+2 
          UU(L,M)=U(L,M,NNN); 
          EUU(L,M)=exp(X(L)*Y(M)*Z(NNN)); 
      end 
  end 
  figure (1); 
  colormap summer; 
  surf(X,Y,UU); 
  xlabel('X values'); 
  ylabel('Y values'); 
  zlabel('UU values'); 
  title('Numerical Solution at z=0.5') 
  figure (2); 
  colormap winter; 
  surf(X,Y,EUU); 
  xlabel('X values'); 
  ylabel('Y values'); 
  zlabel('EUU values'); 
  title('Exact Solution at z=0.5') 