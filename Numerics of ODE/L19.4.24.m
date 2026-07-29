%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Cowell Method
%
% Differential Equation:
%     u'' = w^2*u
%
% Exact Solution:
%     u(t) = cosh(w*t)
%
% Initialization:
%     Exact solution used at t = h
%
% Numerical Method:
%     Two-Step Cowell Method
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


clear all; %clear all variabies in memory
N=10; % No of internal grid points
W=2;
for I=1:8
    H=1.0/N;
    WO=W*W*H*H/12;
disp('Num Soln of IVP using 2-Step Cowels Methods');
disp('     ');
fprintf('N=%4.2f,W=%4.2f, \n',N,W);
disp('     ');
%Formation of Nodal Points
for L=1:N+1
    T(L)=(L-1)*H;
end
%Exact Solution & RHS function
for L=1:N+1
    EU(L)=0.5*(exp(W*T(L))+exp(-W*T(L)));
end
%Initial conditions.
U(1)=EU(1);
U(2)=EU(2);
tic %CPU time starts
%Solution at advanced time point
for L=2:N
    UO=U(L); UUO=U(L-1);
    TO=T(L); TTO=T(L-1);
    U(L+1)=2*U0-UU0+W*W*H*H*U0;
    end
%Calculate Errors
ERU=abs(U(N+1)-EU(N+1));
%Print errors
fprintf('N=%4.2f,ERU=%10.4e\n',N,ERU);
disp('      ');
N=2*N;
end
toc %CPU time ends
%Plot the resulta
figure (1);
plot(T,U,'-g',T,EU,'--r');
xlabel('t-values');
ylabel('Numerical and Exact values u');
title('Exact vs Numerical Solution');
%Multistep (2-step) Numerov Method
%Solution at advanced time point
for L = 2:N
    UO = U(L);
    UUO = U(L - 1);
    TO = T(L);
    TT0= T(L - 1);
    RH= 2*UO -UUO+WO*(UU0+10*U0);
    U(L + 1)=RH/(1 - WO);
end