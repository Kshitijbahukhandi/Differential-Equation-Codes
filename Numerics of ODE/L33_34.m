%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Solution of a System of Initial Value Problems
% using the Third-Order Explicit Runge-Kutta-Nyström Method
%
% Mesh Type:
%     Uniform and Variable Mesh
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
%     Third-Order Explicit Runge-Kutta-Nyström Method
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
ETA=1.2; % Mesh refinement factor 
N=10; 
H=1.0/N; 
disp('Variable Mesh for system of IVPs Using 3rd Nystrom Explicit Method'); 
% Variable Mesh for system of IVPs Using 3rd Nystrom Explicit Method 
disp(' '); 
fprintf('N=%4.2f\n',N); 
N=10.00 
disp(' '); 
for I =1:8 % Loop for mesh refinement 
    SUM=1.0; 
    for L=1:N-1 
SUM=SUM+ETA^L; 
    end 
    H(2)=1.0/SUM; 
    for L=2:N 
H(L+1)=ETA*H(L); 
    end 
    %Formation of grid points 
    T(1)=0.0; 
    for L=2:N+1 
T(L)=T(L-1)+H(L); 
    end 
    T(N+1)=1.0; 
    %Exact solution & RHS function 
    for L=1:N+1 
EU(L)=(exp(-T(L))-exp(-6*T(L)))/5; % Exact solution for u 
EV(L)=(2*exp(-T(L))+3*exp(-6*T(L)))/10; % Exact solution for v 
    end 
    %Initial conditions 
    U(1)=EU(1); 
    V(1)=EV(1); 
    tic % Start CPU time 
    %Solution at advanced time 
    for L=1:N 
H0=H(L+1); % Step size at this iteration 
V0=V(L);  % Previous value of v 
T0=T(L);  % Previous value of t 
U0=U(L);  % Previous value of u 
K11=H0*(-3*U0+2*V0); 
K12=H0*(3*U0-4*V0); 
K21=H0*(-3*(U0+K11)+2*(V0+K12)); % Correction in calculating K21 
K22=H0*(3*(U0+K11)-4*(V0+K12));   
% Correction in calculating K22 
K31=H0*(-3*(U0+K21)+2*(V0+K22)); % Correction in calculating K31 
K32=H0*(3*(U0+K21)-4*(V0+K22));   
% Correction in calculating K32 
U(L+1)=U(L)+(2*K11+3*K21+3*K31)/8; 
V(L+1)=V(L)+(2*K12+3*K22+3*K32)/8; 
end 
%Calculate errors 
ERU=abs(U(N+1)-EU(N+1)); 
ERV=abs(V(N+1)-EV(N+1)); 
%print errors 
fprintf('N=%4.2f, ERU=%10.4e, ERV=%10.4e\n',N,ERU,ERV); 
disp(' '); 
N=2*N; 
end 

toc % End CPU time 
% Elapsed time is 0.037213 seconds. 
%Plot the results 
figure(1); 
plot(T,U,'-g',T,EU,'--r'); 
xlabel('t--values'); 
ylabel('Numerical & Exact values u'); 
title('Exact vs Numerical solution for u'); 
legend('Numerical u', 'Exact u');

figure(2); 
plot(T,V,'-g',T,EV,'--r'); 
xlabel('t--values'); 
ylabel('Numerical and Exact values v'); 
title('Exact vs Numerical solution for v'); 
legend('Numerical v', 'Exact v');