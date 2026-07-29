% Name - KSHITIJ BAHUKHANDI
% Roll Number=SAU/AM(M)/2023/06
% Two-Level Explicit & Implicit Scheme
% Uxx=Ut, % U(x,t)=exp(-pi*pi*t)*sin(pi*x)
% Theta=1/2 (Crandall's Scheme)
% SOUTH ASIAN UNIVERSITY
% Lab File=10
%-----------------------------------------
clear all; % clear all variables in memory
LAM=1.6; X0=0.0; XF=1.0; T0=0.0; TF=1.0;
LL=8;
% LL=16, LL=32, LL=64
H=1.0/LL;
K=LAM*H*H;
%JJ=(TF-T0)/K;
JJ=40; % JJ=160, JJ=640, JJ=2560, %Solution at t=1
disp('Numerical Solution of 1D Heat Conduction Equation');
disp(' ');
fprintf('LL=%4.2f, K=%6.4f, LAM=%6.2f \n', LL, K, LAM);
disp(' ');
% Dimensions
X = zeros(1, LL+1); T = zeros(1, JJ+1);
SUB = zeros(1, LL+1); DIAG = zeros(1, LL+1);
SUP = zeros(1, LL+1); RH = zeros(1, LL+1);
U = zeros(LL+1, JJ+1); EU = zeros(LL+1, JJ+1);
% Formation of Nodal Points
for L = 1:LL+1
X(L) = (L-1)*H;
end
for J = 1:JJ+1
T(J) = (J-1)*K;
end
% Exact Solution & RHS function
for L = 1:LL+1
for J = 1:JJ+1
EU(L, J) = exp(-pi*pi*T(J)) * sin(pi*X(L));
end
end
% Initial conditions
for L = 1:LL+1
U(L, 1) = EU(L, 1);
end
% Boundary conditions
for J = 2:JJ+1
U(1, J) = EU(1, J);
U(LL+1, J) = EU(LL+1, J);
end
% Solution at advanced time level
tic % CPU time starts
for J = 1:JJ
% Solution at Internal grid points at each time level
for L = 2:LL
R0 = U(L, J);
R1 = U(L+1, J);
R2 = U(L-1, J);
SUB(L) = -0.5 * LAM;
DIAG(L) = 1 + LAM;
SUP(L) = -0.5 * LAM;
RH(L) = 0.5 * LAM * R1 + (1 - LAM) * R0 + 0.5 * LAM * R2;
end
RH(2) = RH(2) - SUB(2) * U(1, J+1);
RH(LL) = RH(LL) - SUP(LL) * U(LL+1, J+1);
SUB(2) = 0.0;
SUP(LL) = 0.0;
% Tri-diagonal solver
SUP(2) = SUP(2) / DIAG(2);
RH(2) = RH(2) / DIAG(2);
for L = 3:LL-1
SUP(L) = SUP(L) / (DIAG(L) - SUB(L) * SUP(L-1));
end
for L = 3:LL
RH(L) = (RH(L) - SUB(L) * RH(L-1)) / (DIAG(L) - SUB(L) * SUP(L-1));
end
U(LL, J+1) = RH(LL);
for L = 2:LL-1
U(LL-L+1, J+1) = RH(LL-L+1) - SUP(LL-L+1) * U(LL-L+2, J+1);
end
% Calculate Errors
MAU = 0.0; RMU = 0.0;
for L = 2:LL
ERU = abs(U(L, J+1) - EU(L, J+1));
RMU = RMU + ERU^2;
if MAU <= ERU
MAU = ERU;
end
end
RMU = sqrt(RMU / (LL-1));
% Print errors
fprintf('T=%6.4f, MAU=%10.4e, RMU=%10.4e\n', J*K, MAU, RMU);
% Set New Initial conditions
for L = 1:LL+1
U(L, J) = U(L, J+1);
end
end
disp(' ');
toc % CPU time ends
% Plot the results
figure(1);
colormap spring;
surf(T, X, U);
xlabel('T values');
ylabel('X values');
zlabel('U values');
title('Numerical Solution');
figure(2);
colormap autumn;
surf(T, X, EU);
xlabel('T values');
ylabel('X Values');
zlabel('EU values');
title('Exact Solution');
figure(3);
plot(X, U(:,end), '--g');
title('Numerical Solution at Final Time Step');
figure(4);
plot(X, EU(:,end), '--r');
title('Exact Solution at Final Time Step');
% Plot at t=1
for L = 1:LL+1
US(L) = U(L, JJ+1);
EUS(L) = exp(-pi*pi) * sin(pi*X(L));
end
figure(5);
plot(X, US, '--g', X, EUS, '--r');
xlabel('X values');
ylabel('US values');
title('Numerical Vs Exact Solution at t=1');