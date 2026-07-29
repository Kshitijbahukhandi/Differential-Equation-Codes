%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Cowell Method
%
% Problem Type:
%     Initial Value Problem (IVP)
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
 
 clear all;
 N=10;
 W=2;
 for I=1:8
    H=1.0/N;
    W0=W*W*H/12;
    disp('num sol of ivp using 2-step cowels method');
    disp(' ');
 %formation of nodal points
 for L=1:N+1
        T(L)=(L-1)*H;
 end
 %exact solition of rhs function
 for L=1:N+1
        EU(L)=0.5*(exp(W*T(L))+exp(-W*T(L)));
 end
 %initial conditions
    U(1)=EU(1);
    U(2)=EU(2);
    tic
 %solution at advanced time point
 for L=2:N
        U0=U(L);UU0=U(L-1);
        T0=T(L);TT0=T(L-1);
        U(L+1)=2*U0-UU0+W*W*H*H*U0;
 end
 %calculate errors
    ERU=abs(U(N+1)-EU(N+1));
 %Print errors
    fprintf('N=%4.2f,ERU=%10.4e\n',N,ERU);
    disp(' ')
    N=2*N;
 end
 toc
 %plot the results
 figure(1);
 plot(T,U,'-g',T,EU,'--r');
 xlabel('t-values');
 ylabel('numerical and exact values u');
 title('exact vs numerial solution')