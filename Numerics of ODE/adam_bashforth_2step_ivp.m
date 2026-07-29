%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Adam-Bashforth Method
%
% Numerical Method:
%     Two-Step Adam-Bashforth Explicit Multistep Method
%
% Problem Type:
%     Initial Value Problem (IVP)
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
 for I=1:8
    H=1.0/N;
    disp('Num sol of IVP using 2 step Adam-bashforth Method');
    fprintf('N=4.2f\n',N);
    disp(' ')
 %formation of nodal points
 for L=1:N+1
        T(L)=(L-1)*H;
 end
 %exact solution and rhs function
 for L=1:N+1
        EU(L)=1.0/(1+T(L)*T(L));
 end
 %intial comditions
    U(1)=EU(1);
    U(2)=EU(2);
    U(1)=1.0;
    tic
 %solution at t=h by the help of second order classial rk method
    K1=-2*H*T(1)*U(1)*U(1);
    K2=-2*H*(T(1)+H)*(U(1)+K1)*(U(1)+K1);
    U(2)=U(1)+(K1+K2)/2;
 %Solution at advanced time point
 for L=2:N
        U0=U(L);
        UU0=U(L-1);
        T0=T(L);
        TT0=T(L-1);
        U(L+1)=U0+0.5*H*(-6*T0*U0*U0+2*TT0*UU0*UU0);
 end
 %calulate errors
    ERU=abs(U(N+1)-U(N+1));
 %PRINT ERRORS
    fprintf('N=4.2f,ERU=%10.4e\n',N,ERU);
    disp(' ');
    N=2*N;
 end
 toc
 %plot the results
 figure(1);
 plot(T,U,'g',T,EU,'--r');
 xlabel('t-values');
 ylabel('numerical and exact values u');
 title('exact vs numerial solution')