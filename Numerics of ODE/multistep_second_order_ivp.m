%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a Second-Order Initial Value Problem
% using a Multistep Method
%
% Differential Equation:
%     u'' = 6u^2
%
% Initial Conditions:
%     u(0) = 1
%     u'(0) = -2
%
% Numerical Method:
%     Explicit Two-Step Multistep Method
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
% Clear all variables in memory
N=40;
% No of internal grid points
for I =1:4
H=1.0/N;
disp('Num soln of Non-linear IVP using 2-step Adams Bashward Method');
disp('     ');
fprintf('N=%4.2f, H=%6.4f\n', N, H)
disp('     ');
%Formation of nodal points
for L=1:N+1
    T(L)=(L-1)*H;
end
%Exact solution & RHS function
for L=1:N+1
    EU(L)=1/(1+T(L))^2;
    EV(L)=-2/(1+T(L))^3;
end
%Initial conditions
U(1)=1;
U(2)=EU(2);
V(1)=-2;
V(2)=EV(2);
tic
%Solution at advanced time point
for L=2:N
    U0=U(L);     UU0=U(L-1);
    V0=V(L);     VV0=V(L-1);
    T0=T(L);     TT0=T(L-1);
    %Multi Step Formula
    V(L+1)=(4/3)*V(L)-(1/3)*V(L-1)+(2/3)*H*(6*(U(L))^2);
    U(L+1)=(4/3)*U(L)-(1/3)*U(L-1)+(2/3)*H*(V(L));
    
    
end
%Calculate Errors
ERU=abs(U(L+1)-EU(L+1));
%print errors
fprintf('N=%4.2f, ERU=%10.4e\n',N,ERU);
disp('     ');
N=2*N;
end
toc
%CPU time ends
figure(3);
plot(T,U,'-g',T,EU,'--r');
xlabel('t-values');
ylabel('numerical and exact values u');
title('exact vs numerial solution')
figure (2);
plot(T,V,'-g',T,EV,'--r');
xlabel('t-values');
ylabel('numerical and exact values v');
title('exact vs numerial solution')
