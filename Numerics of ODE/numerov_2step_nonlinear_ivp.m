%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a Nonlinear Initial Value Problem
% using the Two-Step Numerov Method
%
% Numerical Method:
%     Two-Step Numerov Multistep Method
%
% Problem Type:
%     Nonlinear Initial Value Problem (IVP)
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
    disp('num sol of non linear ivp using 2-step numerov method');
    disp(' ');
    fprintf('N=%4.2f,\n',N);
    disp(' ')
 %formation of nodal points
 for L=1:N+1
        T(L)=(L-1)*H;
 end
 %exact solition of rhs function
 for L=1:N+1
        EU(L)=1+exp(2*T(L));
        F(L)=2*exp(2*T(L))-1-exp(4*T(L));
 end
 %initial conditions
    U(1)=EU(1);
    U(2)=EU(2);
    NU(1)=EU(1);
    NU(2)=EU(2);
    tic
 %solution at advanced time point
 for L=2:N
        OLDU(3)=EU(L+1);
        U(3)=OLDU(3);
 for KK=1:600
            F2=F(L);F3=F(L+1);F1=F(L-1);
            F4=H0*(F1+F3+10*F2);
            U1=U(1);U2=U(2);U3=U(3);
            AA=H0*U3*U3-U3+2*U2-U1+H0*(U1*U1+10*U2*U2)+F4;
            DA=2*H0*U3-1;
            U(3)=U(3)-AA/DA;
            NU(L+1)=U(3);
            BIG=abs(U(3)-OLDU(3));
 if BIG<=10^(-12),break;end
 end
        U(1)=U(2);
        U(2)=U(3);
 end
    ERR=abs(NU(N+1)-EU(N+1));
 
%Print errors
    fprintf('N=%4.2f,ERU=%10.4e\n',N,ERR);
    disp(' ')
    N=2*N;
 end
 toc
 %plot the results
 figure(1);
 plot(T,NU,'-g',T,EU,'--r');
 xlabel('t-values');
 ylabel('numerical and exact values u');
 title('exact vs numerial solution')