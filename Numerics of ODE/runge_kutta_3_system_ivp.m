%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Third-Order Explicit Runge-Kutta Method
%
% System of Differential Equations:
%     du/dt = -3u + 2v
%     dv/dt =  3u - 4v
%
% Exact Solutions:
%     u(t) = (exp(-t) - exp(-6t))/5
%     v(t) = (2*exp(-t) + 3*exp(-6t))/10
%
% Numerical Method:
%     Third-Order Explicit Runge-Kutta (RK3) Method
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
ETA=1.2; 
H=1.0/N; 
disp('Variable mesh for system of IVPs- 3rd order classical explicit RK method- variable mesh'); 
% Variable mesh for system of IVPs- 3rd order classical explicit RK method- variable mesh
disp('  '); 
   
fprintf('N=%4.2f\n',N); 
N=10.00 
disp('  '); 
   
for I=1:8 
    %formation of variable mesh 
    SUM=1.0; 
    for L=1:N-1 
        SUM=SUM+ETA^L; 
    end 
    H(2)=1.0/SUM; 
    for L=2:N 
        H(L+1)=ETA*H(L); 
    end 
    %formation of grid points 
    T(1)=0.0; 
    for L=2:N+1 
        T(L)=T(L-1)+H(L); 
    end 
    T(N+1)=1.0; 
    %exact solution and rhs function 
    for L=1:N+1 
        EU(L)=(exp(-T(L))-exp(-6*T(L)))/5; 
        EV(L)=(2*exp(-T(L))+3*exp(-6*T(L)))/10; 
    end 
    %initial conditions 
    %U(1)=EU(1); 
    %V(1)=EV(1); 
    U(1)=0.0; 
    V(1)=0.5; 
    tic 
    %solution at advanced time 
    for L=1:N 
        T0=T(L); 
        U0=U(L); 
        V0=V(L); 
        H0=H(L+1); 
        K11=H0*(-3*U0+2*V0); 
        K12=H0*(3*U0-4*V0); 
        K21=H0*(-3*(U0+2*K11/2)+2*(V0+2*K12/2)); 
        K22=H0*(3*(U0+2*K11/2)-4*(V0+2*K12/2)); 
        K31=H0*(-3*(U0-K11+2*K21)+2*(V0-K12+2*K22)); 
        K32=H0*(3*(U0-K11+2*K21)-4*(V0-K12+2*K22)); 
        U(L+1)=U(L)+(K11+4*K21+K31)/6; 
        V(L+1)=V(L)+(K12+4*K22+K32)/6; 
    end 
    %calculate errors 
    ERU=abs(U(N+1)-EU(N+1)); 
    ERV=abs(V(N+1)-EV(N+1)); 
    %print errors 
    fprintf('N=%4.2f,ERU=%10.4e\n',N,ERU,ERV); 
    disp('  '); 
    N=2*N; 
end 
toc 
% Elapsed time is 0.021920 seconds. 
%plot the results 
figure(1); 
plot(T,U,'--g'); 
xlabel('t-values'); 
ylabel('Numerical values u'); 
title('numerical solution') 

figure(2); 
plot(T,EU,'--r'); 
xlabel('t-values'); 
ylabel('exact values u'); 
title('exact solution')

figure(3); 
plot(T,U,'-g',T,EU,'--r'); 
xlabel('t-values'); 
ylabel('numerical and exact values u'); 
title('exact vs numerical solution') 