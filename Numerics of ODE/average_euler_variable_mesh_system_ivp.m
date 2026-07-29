%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Average Euler Method on a Variable Mesh
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
%     Average Euler Method with Variable Mesh
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
ETA=1.2 
ETA = 1.2000 
N=10; 
for I=1:8 
    disp('variable mesh method for non-linear ivp - average euler method'); 
    disp(' '); 
    fprintf('N=%4.2f\n',N); 
    disp(' '); 
    %formation of variable mesh lengths 
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
        EU(L)=1.0/(1+T(L)*T(L)); 
    end 
    %initial conditions 
    U(1)=EU(1); 
    U(2)=U(1)-H(2)*H(2); 
    tic %CPU time starts 
    %solution at advanced time point 
    for L=2:N 
        U(L+1)=U(L-1)-2*(H(L)+H(L+1))*T(L)*U(L)^2; 
    end 
    N=2*N; 
end

toc 
% Elapsed time is 0.027950 seconds. 
%plot the results 
figure(1); 
plot(T,U,'--g'); 
xlabel('t-values'); 
ylabel('numerical values u'); 
title('numerical solution') 

figure(2); 
plot(T,EU,'--r'); 
xlabel('t-values'); 
ylabel('exact values u'); 
title('exact solution');

figure(3); 
plot(T,U,'-g',T,EU,'--r'); 
xlabel('t-values'); 
ylabel('numerical and exact values u'); 
title('exact vs numerical solution'); 

figure(4); 
plot(T,U,'b+');hold on 
plot(T,EU,'-r');hold off 
xlabel('t-values'); 
ylabel('numerical and exact values u'); 
title('exact vs numerical solution'); 