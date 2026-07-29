%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second-Order Numerical Solution of a Linear Two-Point Boundary Value Problem
%
% Problem Type:
%     One-Dimensional Linear Elliptic Boundary Value Problem
%
% Differential Equation:
%     u'' = u - 1
%
% Boundary Conditions:
%     u(0) - u'(0) = 1
%     u(1) + u'(1) = 1 + 2e
%
% Boundary Discretization:
%     Symmetric Discretization of Mixed Boundary Conditions
%
% Exact Solution:
%     u(x) = 1 + exp(x)
%
% Numerical Method:
%     Second-Order Finite Difference Method
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
LL=80;
H=1.0/LL;
disp('Num Soln for 1D linear Elliptic Equation with mixed BCs');
disp('    ');
fprintf('LL=%4.2f\n',LL);
disp('   ')
%Formation of Nodal Points
for L=1:LL+1
    X(L)=(L-1)*H;
end
% Exact Solution & RHS function
for L=1:LL+1
    EU(L)=1+exp(X(L));
    F(L)=-1.0;
end
% Boundary condtions
%U(1)=EU(1);
%U(LL+1)=EU(LL+1);
tic % cpu time starts
DIAG(1)=-(2+2*H+H*H);
SUP(1)=2.0;
RH(1)=-2*H-H*H;
for L=2:LL
    F0=F(L); 
    SUB(L)=1.0;
    DIAG(L)=-(2+H*H);
    SUP(L)=1.0;
    RH(L)=H*H*F0;
end
SUB(LL+1)=2.0;
DIAG(LL+1)=-(2+2*H+H*H);
RH(LL+1)=-2*H*(1+2*exp(1))-H*H;
SUB(1)=0.0;
SUP(LL+1)=0.0;
% Tri-diagonal solver 
SUP(1)=SUP(1)/DIAG(1);
RH(1)=RH(1)/DIAG(1);
for L=2:LL
    SUP(L)=SUP(L)/(DIAG(L)-SUB(L)*SUP(L-1));
end
for L=2:LL+1
    RH(L)=(RH(L)-SUB(L)*RH(L-1))/(DIAG(L)-SUB(L)*SUP(L-1));
end
U(LL+1)=RH(LL+1);
for L=1:LL
    U(LL-L+1)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+1+1);
end
% Calculate Errors
MAU=0.0;
RMU=0.0;
for L=1:LL+1
    ERU=abs(U(L)-EU(L));
    RMU=RMU+ERU^2;
    if MAU<=ERU, MAU=ERU;  end
end
RMU=sqrt(RMU/(LL+1));
% Print errors 
fprintf('MAU=%10.4e, RMU=%10.4e\n',MAU,RMU);
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