%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Numerov Method
%
% Differential Equation:
%     u'' = u^2 + f(t)
%
% Exact Solution:
%     u(t) = 1 + exp(2t)
%
% Initialization:
%     Exact solution used at t = h
%
% Numerical Method:
%     Two-Step Numerov Method
%
% Language:
%     MATLAB
%
% Developed as part of the coursework for
% Numerical Methods for Ordinary Differential Equations
% M.Sc. Applied Mathematics
% South Asian University
%
% Author: Kshitij Bahukhandi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
clear all; 
N=10; 
for I=1:7 
    H=1.0/N; 
    H0=H*H/12; 
    disp('Num Soln of Nonlinear IVP using 2-step Numerov Method'); 
    disp('  '); 
    fprintf('N=%4.2f,\n',N); 
    disp('  '); 
    % Formation of Nodal Points 
    for L=1:N+1 
        T(L)=(L-1)*H; 
    end 
    % Exact Solution & RHS Function  
    for L=1:N+1 
        EU(L)=1+exp(2*T(L)); 
        F(L)=2*exp(2*T(L))-1-exp(4*T(L)); 
    end 
    % Initial Conditions 
    U(1)=EU(1); 
    U(2)=EU(2); 
    NU(1)=EU(1); 
    NU(2)=EU(2); 
     
    tic % CPU time starts 
     
    % Solution at advanced time point 
    for L=2:N 
        OLDU(3)=EU(L+1); 
        U(3)=OLDU(3); 
        for KK=1:600 
            F2=F(L); F3=F(L+1); F1=F(L-1); 
            F4=H0*(F1+F3+10*F2); 
            U1=U(1); U2=U(2); U3=U(3); 
            AA=H0*U3*U3-U3+2*U2-U1+H0*(U1*U1+10*U2*U2)+F4; 
            DA=2*H0*U3-1; 
            U(3)=U(3)-AA/DA; 
            NU(L+1)=U(3); 
            BIG=abs(U(3)-OLDU(3)); 
            if BIG <=10^(-12), break; end 
            OLDU(3)=U(3); 
        end 
        U(1)=U(2); 
        U(2)=U(3); 
    end 
     
   ERR=abs(NU(N+1)-EU(N+1)); 
    
   % Print Errors 
    fprintf('N=%4.2f, ERR=%10.4e\n',N,ERR); 
    disp('  '); 
    N=2*N; 
end 
  
toc % CPU time ends 

% Plot the results 
figure(1); 
plot(T,NU,'-g',T,EU,'--r'); 
xlabel('t-values'); 
ylabel('Numerical and Exact values u'); 
title('Exact vs Numerical Solution') 