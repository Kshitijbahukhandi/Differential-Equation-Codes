%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth-Order Numerical Solution of a Linear Two-Point Boundary Value Problem
%
% Problem Type:
%     One-Dimensional Linear Elliptic Boundary Value Problem
%
% Differential Equation:
%     u'' - ?²u = f(x)
%
% Exact Solution:
%     u(x) = exp(x^2)
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
H=1.0/LL;
LAM=10;
H0=LAM*LAM*H*H/12;
disp('Fourth order Num Soln for 1D linear Elliptic Equation');
disp('    ');
fprintf('LL=%4.2f,LAM=%8.4f\n',LL,LAM);
disp('   ')
%Formation of Nodal Points
for L=1:LL+1
    X(L)=(L-1)*H;
end
% Exact Solution & RHS function
for L=1:LL+1
    EU(L)=exp(X(L)^2);
    F(L)=(2-LAM*LAM+4*X(L)*X(L))*EU(L);
end
% Boundary condtions
U(1)=EU(1);
U(LL+1)=EU(LL+1);
tic % cpu time starts
for L=2:LL
    F0=F(L); F1=F(L+1); F2=F(L-1);
    F3=H*H*(F1+F2+10*F0)/12;
    U0=U(L);
    U1=U(L+1);
    U2=U(L-1);
    SUB(L)=1-H0;
    DIAG(L)=-(2+10*H0);
    SUP(L)=1-H0;
    RH(L)=F3;
end
RH(2)=RH(2)-SUB(2)*U(1);
RH(LL)=RH(LL)-SUP(LL)*U(LL+1);
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
U(LL)=RH(LL);
for L=2:LL-1
    U(LL-L+1)=RH(LL-L+1)*U(LL-L+2);
end
U(LL)=RH(LL);
for L=2:LL-1
    U(LL-L+1)=RH(LL-L+1)-SUP(LL-L+1)*U(LL-L+2);
end
% Calculate Errors
MAU=0.0;
RMU=0.0;
for L=2:LL
    ERU=abs(U(L)-EU(L));
    RMU=RMU+ERU^2;
    if MAU<=ERU, MAU=ERU;  end
end
RMU=sqrt(RMU/(LL-1));
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