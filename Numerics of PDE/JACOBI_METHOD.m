%-------------------------------------------------------------------------- 
% NAME:- KSHITIJ BAHUKHANDI
% ENROLLMENT NUMBER:- SAU\AM(M)\2023\06
%MSC IN APPLIED MATHEMATICS (SEMESTER 3).
%LAB = 01.    %MATLAB FILE:- 01.      %DATE=06.08.2024
%JACOBI METHOD.
%NUMERICAL SOLUTION OF 2D POISSON EQUATION.
%Uxx+Uyy=f(x,y), %U(x,y)=exp(x*y).
%SOUTH ASIAN UNIVERSITY.
%--------------------------------------------------------------------------
clear all;  
LL=8; 
MM=LL; 
H=1.0/(LL+1); 
disp('Numerical Solution of 2D Poissons Equation-Jacobi Method'); 
disp('   '); 
fprintf('LL=%4.2f\n',LL); 
disp('   '); 
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
% Boundary Conditions  
for L=1:LL+2 
for M=1:MM+2 
U(L,1)=EU(L,1); 
U(L,MM+2)=EU(L,MM+2); 
U(1,M)=EU(1,M); 
U(LL+2,M)=EU(LL+2,M); 
end 
end 
% Initial and First Approximations 
for L=1:LL+2 
for M=1:MM+2 
OLDU(L,M)=EU(L,M); 
U(L,M)=OLDU(L,M); 
end 
end 
tic % CPU time starts
% Iteration Starts  
 
for KK=1:6000 
    % For J-iteration put U0=OLDU(L,M),... 
    %--------------------------------- 
    for M=2:MM+1 
        for L=2:LL+1 
            F0=F(L,M); 
            U0=OLDU(L,M); 
            U1=OLDU(L+1,M); 
            U2=OLDU(L-1,M); 
            U3=OLDU(L,M+1); 
            U4=OLDU(L,M-1); 
            U(L,M)=(U1+U2+U3+U4-H*H*F0)/4; 
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
     
    fprintf('KK=%4.0f, MAU=%10.4e, RMU=%10.4e\n',KK,MAU,RMU); 
    % Error Tolerance  
    BIG=0.0; 
    for L=2:LL+1 
        for M=2:MM+1 
             
            FU=abs(U(L,M)-OLDU(L,M)); 
            if FU>=BIG, BIG=FU; end 
        end 
    end 
     
    if BIG<=10^(-12),break; end 
    for L=2:LL+1 
        for M=2:MM+1 
            OLDU(L,M)=U(L,M); 
        end 
    end 
    end 
disp('  '); 
toc % CPU time ends 
colormap spring 
figure(1); 
surf(X,Y,EU); 
xlabel('X values'); 
ylabel('Y values'); 
zlabel('U values'); 
title('numerical solution'); 
figure(2); 
surf(X,Y,EU); 
title('Exact Solution'); 
xlabel('X values'); 
ylabel('Y values'); 
zlabel('EU values'); 