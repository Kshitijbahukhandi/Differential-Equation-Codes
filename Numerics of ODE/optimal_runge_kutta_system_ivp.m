%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Optimal Runge-Kutta Method
%
% System of Differential Equations:
%     du/dt = -3u + 2v
%     dv/dt =  3u - 4v
%
% Exact Solutions:
%     u(t) = (exp(-t)-exp(-6t))/5
%     v(t) = (2*exp(-t)+3*exp(-6t))/10
%
% Numerical Method:
%     Optimal Runge-Kutta Method
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
disp('Variable mesh for system of IVPs- Taylor series method- variable mesh'); 
% Variable mesh for system of IVPs- Taylor series method- variable mesh 
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
        DU1=-3*U0+2*V0; 
        DV1=3*U0-4*V0; 
        DU2=-3*DU1+2*DV1; 
        DV2=3*DU1-4*DV1; 
        DU3=-3*DU2+2*DV2; 
        DV3=3*DU2-4*DV2; 
        DU4=-3*DU3+2*DV3; 
        DV4=3*DU3-4*DV3; 
        %Second order approximation 
        U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2; 
        V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2; 
        %Third order approxiamtion  
        U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*DU3/6; 
        V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2+H0*H0*H0*DV3/6; 
        %fourth order approximation 
        
U(L+1)=U(L)+H0*DU1+0.5*H0*H0*DU2+H0*H0*H0*(DU3/6)+H0*H0*H0*H0*(DU4/24); 
        
V(L+1)=V(L)+H0*DV1+0.5*H0*H0*DV2+H0*H0*H0*(DV3/6)+H0*H0*H0*H0*(DV4/24); 
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
% Elapsed time is 0.032059 seconds. 
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