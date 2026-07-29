%----------------------------------------------------- 
% NAME:- KSHITIJ BAHUKHANDI
% South Asian University 
% ENROLLMENT NUMBER:- SAU\AM(M)\2023\06
% 3rd Sem, Numerics of PDE 
% Lab File 04 
% Numerical solution of 2D Nonlinear Elliptic Equation using o(h^2) method 
% Uxx+Uyy=BETA*U*U+f(x,y), %U(x,y)=exp(x)*sin(pi*y) 
%----------------------------------------------------- 
clear all; 
BETA=10.0; 
LL=7; 
MM=LL; 
H=1.0/(LL+1); 
disp('Numerical Solution of 2D Nonlinear Elliptic Equation'); 
disp('  '); 
fprintf('LL=%4.2f,  BETA=%6.0f\n',LL,BETA); 
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
EU(L,M)=exp(X(L))*sin(pi*Y(M)); 
F(L,M)=(1-pi*pi)*EU(L,M)-BETA*EU(L,M)*EU(L,M); 
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
OLDU(L,M)=EU(L,M); 
U(L,M)=OLDU(L,M); 
end 
end 
tic  
  
for KK=1:6000 
     
    %----------------------------------------------------------------------
     
    for M=2:MM+1 
        for L=2:LL+1 
            F0=F(L,M); 
            FF=H*H*F0; 
            U0=U(L,M); 
            U1=U(L+1,M); 
            U2=U(L-1,M); 
            U3=U(L,M+1); 
            U4=U(L,M-1); 
             
            AA=BETA*H*H*U0*U0+4*U0+FF-U1-U2-U3-U4; 
            DA=2*BETA*H*H*U0+4; 
            U(L,M)=U(L,M)-AA/DA; 
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
     
    % Print Errors 
    fprintf('KK=%4.0f,  MAU=%10.4e  , RMU=%10.4e\n',KK,MAU,RMU); 
     
     
    % Error Tolerance  
    BIG=0.0; 
    for L=2:LL+1 
        for M=2:MM+1 
            FU=abs(U(L,M)-OLDU(L,M)); 
            if FU>=BIG, BIG=FU; end 
        end 
    end 
     
    if BIG<=10^(-12), break; end 
    for L=2:LL+1 
        for M=2:MM+1 
OLDU(L,M)=U(L,M); 
end 
end 
end 
disp('   '); 
toc 
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