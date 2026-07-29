%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth-Order Numerical Solution of a Nonlinear Two-Point Boundary Value Problem
%
% Problem Type:
%     One-Dimensional Nonlinear Elliptic Boundary Value Problem
%
% Differential Equation:
%     u'' = 6u^2
%
% Boundary Conditions:
%     u(0) - u'(0) = 3
%     u(1) + u'(1) = 0
%
% Exact Solution:
%     u(x) = 1/(1 + x)^2
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
LL=20;
H=1.0/(LL+1);
H0=H*H/12; H1=1+H; H2=1+0.5*H; H3=3*H*H/4; H4=3*H/2;
disp('Fourth order Num Soln for 1D Non-linear Elliptic BVP with mixed BCs');
disp('    ');
fprintf('LL=%4.2f\n',LL);
disp('   ');
%Formation of Nodal Points
for L=1:LL+2
    X(L)=(L-1)*H;
end
% Exact Solution & RHS function
for L=1:LL+2
    EU(L)=1/((1+X(L))^2);
end
%Initial and First Approximations
for L=1:LL+2
    OLDU(L)=EU(L);
    U(L)=OLDU(L);
end
tic % cpu time starts
%Iteration starts
for KK=1:6000
    %--------------
    A0=H2*U(1)+H3*U(1)*U(1)-H4;
    AA=H1*U(1)+H*H*U(1)*U(1)+2*H*H*A0*A0-3*H-U(2);
    DA=H1+2*H*H*U(1)+4*H*H*A0*(H2+2*H3*U(1));
    U(1)=U(1)-AA/DA;
for L=2:LL+1
    U0=U(L);
    U1=U(L+1);
    U2=U(L-1);
    BB=10*H0*U0+2*U0+H0*(U1*U1+U2*U2)-U1-U2;
    DB=20*H0*U0+2;    
    U(L)=U(L)-BB/DB;
end
C0=H2*U(LL+2)+H3*U(LL+2)*U(LL+2);
CC=H1*U(LL+2)+H*H*U(LL+2)+2*H*H*C0*C0-U(LL+2);
DC=H1+2*H*H*U(LL+2)+4*H*H*C0*(H2+2*H3*U(LL+2));
U(LL+2)=U(LL+2)-CC/DC;
end
% Calculate Errors
MAU=0.0;
RMU=0.0;
for L=1:LL+2
    ERU=abs(U(L)-EU(L));
    RMU=RMU+ERU^2;
    if MAU<=ERU, MAU=ERU; end
end
RMU=sqrt(RMU/LL+2);
% Print errors 
fprintf('KK=%4.0f, MAU=%10.4e, RMU=%10.4e\n',KK,MAU,RMU);
% Error Tolerance
BIG=0.0;
for L=1:LL+2
    FU=abs(U(L)-OLDU(L));
    if FU >= BIG, BIG=FU; end
    if BIG <= 10^(-12), break; end
    for L=1:LL+2
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