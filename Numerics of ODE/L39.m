%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Adams-Moulton Method
%
% Differential Equation:
%     du/dt = -2*t*u^2
%
% Exact Solution:
%     u(t) = 1/(1+t^2)
%
% Initialization:
%     Third-Order Classical Explicit Runge-Kutta (RK3) Method
%
% Numerical Method:
%     Two-Step Adams-Moulton Method
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
% Initial value for N 
N=10;  
% Loop for various values of N 
for I=1:8 
    H=1.0/N; 
    disp('Num Soln of Non-Linear IVP using 2-step Adams-Moulton Method'); 
    disp('  '); 
     
    % Formation of Nodal Points 
    for L=1:N+1 
        T(L)=(L-1)*H; 
    end 
     
    % Exact Solution & RHS function  
    for L=1:N+1 
        EU(L)=1.0/(1+T(L)*T(L)); 
    end 
     
    % Initial Conditions 
    U(1)=EU(1); 
    NU(1)=EU(1); 
    tic % CPU time starts 
     
    % Solution at t=h using 3rd order classical explicit RK-Method  
    K1=-2*H*T(1)*U(1)*U(1); 
    K2=-2*H*(T(1)+H/2)*(U(1)+K1/2)*(U(1)+K1/2); 
    K3=-2*H*(T(1)+H)*(U(1)-K1+2*K2)*(U(1)-K1+2*K2); 
    U(2)=U(1)+(K1+4*K2+K3)/6; 
    NU(2)=EU(2); 
     
    % Solution at advanced time point  
    for L=2:N 
        OLDU(3)=EU(L+1); 
        U(3)=OLDU(3); 
        for KK=1:600 
            T1=T(L-1); T2=T(L); T3=T(L+1); 
            U1=U(1); U2=U(2); U3=U(3); 
            AA=(5*H/6)*T3*U3*U3+U3-U2+H*(8*T2*U2*U2-T1*U1*U1)/6; 
            DA=(5*H/3)*T3*U3+1; 
            U(3)=U(3)-AA/DA; 
NU(L+1)=U(3); 
BIG=abs(U(3)-OLDU(3)); 
if BIG<=10^(-12), break; 
end 
OLDU(3)=U(3); 
end 
% Set new Initial Conditions 
U(1)=U(2); 
U(2)=U(3); 
end 
% Calculate error 
ERR=abs(NU(N+1)-EU(N+1)); 
% Print errors 
fprintf('N=%4.2f, ERR=%10.4e/n',N,ERR); 
disp('  '); 
% Update N 
N=2*N; 
end 

%%
toc % CPU Time ends 
% Elapsed time is 0.038832 seconds. 
% Plot the results 
figure(1); 
plot(T,NU,'-g',T,EU,'--r'); 
xlabel('t-values'); 
ylabel('Numerical and Exact Values u'); 
title('Exact vs Numerical Solution');
%