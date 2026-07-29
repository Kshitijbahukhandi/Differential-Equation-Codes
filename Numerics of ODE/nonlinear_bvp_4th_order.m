%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth-Order Numerical Solution of a Nonlinear Two-Point Boundary Value Problem
%
% Problem Type:
%     One-Dimensional Nonlinear Elliptic Boundary Value Problem
%
% Differential Equation:
%     u'' = ?u² + f(x)
%
% Exact Solution:
%     u(x) = sin(?x)
%
% Numerical Method:
%     Fourth-Order Finite Difference Method
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

clear all; % clear all variables in memory
LL=40;
H=1.0/(LL+1);
alpha=10.0;
H0=H*H/12;
disp('Fourth order Num Soln for 1D Non-linear Elliptic Equation');
disp('    ');
fprintf('LL=%4.2f,alpha=%8.4f\n',LL,alpha);
disp('   ')
%Formation of Nodal Points
for L=1:LL+2
    X(L)=(L-1)*H;
end
% Exact Solution & RHS function
for L=1:LL+2
    EU(L)=sin(pi*X(L));
    F(L)=-pi*pi*EU(L)-alpha*EU(L)*EU(L);
end
% Boundary condtions
U(1)=EU(1);
U(LL+2)=EU(LL+2);
%Initial and First Approximations
for L=2:LL+1
    OLDU(L)=EU(L);
    U(L)=OLDU(L);
end
tic % cpu time starts
%Iteration starts
for KK=1:600
    %--------------
for L=2:LL+1
    F0=F(L); F1=F(L+1); F2=F(L-1);
    F3=F1+F2+10*F0;
    U0=U(L);
    U1=U(L+1);
    U2=U(L-1);
    A0=10*alpha*H0;
    C0=H0*(alpha*(U1*U1+U2*U2)+F3)-U1-U2;
    AA=A0*U0*U0+2*U0+C0;
    DA=2*A0*U0+2;
    U(L)=U(L)-AA/DA;
end
end
% Calculate Errors
MAU=0.0;
RMU=0.0;
for L=2:LL+1
    ERU=abs(U(L)-EU(L));
    RMU=RMU+ERU^2;
    if MAU<=ERU, MAU=ERU; end
end
RMU=sqrt(RMU/LL);
% Print errors 
fprintf('KK=%4.0f, MAU=%10.4e, RMU=%10.4e\n',KK,MAU,RMU);
% Error Tolerance
BIG=0.0;
for L=2:LL+1
    FU=abs(U(L)-OLDU(L));
    if FU >= BIG, BIG=FU; end
    if BIG <= 10^(-12), break; end
    for L=2:LL+1
        OLDU(L)=U(L);
    end
end    
disp('  ');
toc  % CPU time ends
% Plot the results 
figure(1);
plot(X,U,'--g');
xlabel('x-values');
ylabel('Numerical values u');
title('Numerical Solution')
figure(2);
plot(X,EU,'--r');
xlabel('x-values');
ylabel('Exact values u');
title('Exact Solution')