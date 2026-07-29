%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of an Initial Value Problem
% using the Two-Step Milne-Simpson Method
%
% Numerical Method:
%     Two-Step Milne-Simpson Multistep Method
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
    disp('Num sol of non-linear IVP using 2 step Adam-Moulton Method');
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
    NU(1)=EU(1);
    tic
 %solution at t=h using 3rd order classical rk method
    K1=-2*H*T(1)*U(1)*U(1);
    K2=-2*H*(T(1)+H/2)*(U(1)+K1/2)*(U(1)+K1/2);
    K3=-2*H*(T(1)+H/2)*(U(1)+K2/2)*(U(1)+K2/2);
    K4=-2*H*(T(1)+H)*(U(1)+K3)*(U(1)+K3);
    U(2)=U(1)+(K1+2*K2+2*K3+K4)/6;
    NU(2)=EU(2);
 %Solution at advanced time point
 for L=2:N
        OLDU(3)=EU(L+1);
        U(3)=OLDU(3);
 for KK=1:600
            T1=T(L-1);
            T2=T(L);
            T3=T(L+1);
            U1=U(1);
            U2=U(2);
            U3=U(3);
            AA=(2*H/3)*T3*U3*U3+U3-U1+(2*H/3)*(4*T2*U2*U2+T1*U1*U1);
            DA=(4*H/3)*T3*U3+1;
            U(3)=U(3)-AA/DA;
            NU(L+1)=U(3);
            BIG=abs(U(3)-OLDU(3));
 if BIG<=10^(-12),break;end
             OLDU(3)=U(3);
 end
 %Set new initial cinditions
        U(1)=U(2);
        U(2)=U(3);
 end
    ERR=abs(NU(N+1)-EU(N+1));
 %PRINT ERRORS
    fprintf('N=4.2f,ERR=%10.4e\n',N,ERR);
    disp(' ');
    N=2*N;
 end
 toc
 %plot the results
 figure(1);
 plot(T,NU,'g',T,EU,'--r');
 xlabel('t-values');
 ylabel('numerical and exact values u');
 title('exact vs numerial solution')