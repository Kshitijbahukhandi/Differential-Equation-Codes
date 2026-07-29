%-------------------------------- 
% NAME:- KSHITIJ BAHUKHANDI
% South Asian University 
% ENROLLMENT NUMBER:- SAU\AM(M)\2023\06
% 3rd Sem, Numerics of PDE 
% Lab File 06 
% Numerical Solution of 2D Poissons Equation using o(h^4) method 
% Uxx+Uyy=f(x,y), % U(x,y)=exp(x*y) 
%----------------------------- 
clear all;  % clear all variables in memory 
LL=7; 
MM=LL; 
H=1.0/(LL+1); 
disp('Numerical Solution of 2D Poissons Equation'); 
disp('    '); 
fprintf('LL=%4.2f \n',LL); 
disp('  '); 
% Formation of Nodal Points 
for L=1:LL+2 
for M=1:MM+2 
X(L)=(L-1)*H; 
Y(M)=(M-1)*H; 
end 
end 
% Exact Solution & RHS function 
for L=1:LL+2 
for M=1:MM+2 
EU(L,M)=exp(X(L)*Y(M)); 
F(L,M)=(X(L)^2+Y(M)^2)*EU(L,M); 
end 
end 
% Boundary conditions 
for L=1:LL+2 
for M=1:MM+2 
U(L,1)=EU(L,1); 
U(L,MM+2)=EU(L,MM+2); 
U(1,M)=EU(1,M); 
U(LL+2,M)=EU(LL+2,M); 
end 
end 
% Initial and First Approximations 
for L=2:LL+1 
for M=2:MM+1 
OLDU(L,M)=0.0; 
U(L,M)=OLDU(L,M); 
end 
end 
tic  
% CPU time starts 
% Iteration starts  
for KK=1:6000 
%----------------  
        for M=2:MM+1 
            for L=2:LL+1 
                F0=F(L,M); 
                F1=F(L+1,M); 
                F2=F(L-1,M); 
                F3=F(L,M+1); 
                F4=F(L,M-1); 
                FF=0.5*H*H*(F1+F2+F3+F4+8*F0); 
                U0=U(L,M); 
                U1=U(L+1,M); 
                U2=U(L-1,M); 
                U3=U(L,M+1); 
                U4=U(L,M-1); 
                U5=U(L+1,M+1); 
                U6=U(L+1,M-1); 
                U7=U(L-1,M+1); 
                U8=U(L-1,M-1); 
                U(L,M)=(U5+U6+U7+U8+4*(U1+U2+U3+U4)-FF)/20; 
            end 
        end 
        % Calculate Errors 
        MAU=0.0; 
        RMU=0.0; 
        for L=2:LL+1 
            for M=2:MM+1 
                ERU=abs(U(L,M)-EU(L,M)); 
                RMU=RMU+ERU^2; 
                if MAU<=ERU, MAU=ERU; end 
            end 
        end 
                RMU=sqrt(RMU/(LL^2)); 
                % MAU=max(max(abs(U-EU))); 
     % Print errors 
      fprintf('KK=%4.0f,  MAU=%10.4e, RMU=%10.4e\n',KK,MAU,RMU); 
      % Error Tolerance 
      %--------------------------- 
      ERR=0.0; 
      for L=2:LL+1 
            for M=2:MM+1 
        FU=abs(U(L,M)-OLDU(L,M)); 
            if FU >= ERR, ERR=FU; end 
            end 
      end   
            if ERR <= 10^(-12), break; end 
          %    for L=2:LL+1 
           %  for M=2:MM+1  
         % OLDU(L,M)=U(L,M); 
         %    end 
         %     end  
         OLDU=U; 
    end 
 disp('  '); 
toc  
% CPU time ends 
figure(1); 
surf(X,Y,U); 
xlabel('X values'); 
ylabel('Y values'); 
zlabel('U values'); 
title('Numerical Solution'); 
figure(2); 
surf(X,Y,EU); 
title('Exact Solution');         
xlabel('X values'); 
ylabel('Y values'); 
zlabel('EU values'); 