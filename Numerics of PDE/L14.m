% Numerical Solution of 2D Nonlinear Parabolic Equation
% Uxx + Uyy = Ut + alpha * U * U + f(x,y,t)
% U(x,y,t) = exp(-t) * cos(pi*x) * cos(pi*y)
% Name - KSHITIJ BAHUKHANDI
% Enrollment Number: SAU/AM(M)/2023/06
% SOUTH ASIAN UNIVERSITY
% Lab File: 14
%-----------------------------
clear all; % Clear all variables in memory
% Parameters
alpha = 10.0;
LAM = 0.8;
LL = 8;
MM = LL; % MM = LL;
JJ = 80; % JJ = 80; % Adjust as needed for more accuracy
H = 1.0 / LL;
K = LAM * H^2;
% Display the problem statement
disp('Numerical Solution of 2D Nonlinear Parabolic Equation');
disp(' ');
fprintf('LL=%4.2f, K=%6.4f, LAM=%6.2f \n', LL, K, LAM);
disp(' ');
% Dimensions
X = zeros(1, LL+1);
Y = zeros(1, MM+1);
T = zeros(1, JJ+1);
TT = zeros(1, JJ+1);
U = zeros(LL+1, MM+1, JJ+1);
EU = zeros(LL+1, MM+1, JJ+1);
OLDU = zeros(LL+1, MM+1, JJ+1);
% Formation of Nodal Points
for L = 1:LL+1
    for M = 1:MM+1
for J = 1:JJ+1
X(L) = (L-1) * H;
Y(M) = (M-1) * H;
T(J) = (J-1) * K;
TT(J) = T(J) + 0.5 * K;
end
    end
end
% Exact Solution & RHS function
for L = 1:LL+1
for M = 1:MM+1
for J = 1:JJ+1
EU(L, M, J) = exp(-T(J)) * cos(pi * X(L)) * cos(pi * Y(M));
ETT(L, M, J) = exp(-TT(J)) * cos(pi * X(L)) * cos(pi * Y(M));
F(L, M, J) = (1 - 2 * pi^2) * ETT(L, M, J) - alpha * ETT(L, M, J)^2;
end
end
end
% Initial conditions
for L = 1:LL+1
for M = 1:MM+1
U(L, M, 1) = EU(L, M, 1);
end
end
% Boundary conditions
for J = 1:JJ
for L = 1:LL+1
for M = 1:MM+1
U(1, M, J+1) = EU(1, M, J+1);
U(LL+1, M, J+1) = EU(LL+1, M, J+1);
U(L, 1, J+1) = EU(L, 1, J+1);
U(L, MM+1, J+1) = EU(L, MM+1, J+1);
end
end
end
% Solution at advanced time level
tic % CPU time starts
for J = 1:JJ
% Initial & First Approx
for L = 2:LL
for M = 2:MM
OLDU(L, M, J+1) = EU(L, M, J+1);
U(L, M, J+1) = OLDU(L, M, J+1);
end
end
for KK = 1:600 % Iteration Starts
% Solution at Internal grid points at each time level
for L = 2:LL
for M = 2:MM
F0 = F(L, M, J);
R0 = U(L, M, J);
R1 = U(L+1, M, J);
R2 = U(L-1, M, J);
R3 = U(L, M+1, J);
R4 = U(L, M-1, J);
R5 = R1 + R2 - 2 * R0;
R6 = R3 + R4 - 2 * R0;
U0 = U(L, M, J+1);
U1 = U(L+1, M, J+1);
U2 = U(L-1, M, J+1);
U3 = U(L, M+1, J+1);
U4 = U(L, M-1, J+1);
U5 = U1 + U2;
U6 = U3 + U4;
P0 = 0.25 * alpha * K;
P1 = 1 + 2 * LAM + 2 * P0 * R0;
P2 = K * F0 + P0 * R0^2 - R0 - 0.5 * LAM * (R5 + R6 + U5 + U6);
AA = P0 * U0^2 + P1 * U0 + P2;
DA = 2 * P0 * U0 + P1;
U(L, M, J+1) = U(L, M, J+1) - AA / DA;
end
end
% Calculate Errors
MAU = 0.0;
RMU = 0.0;
for L = 2:LL
for M = 2:MM
ERU = abs(U(L, M, J+1) - EU(L, M, J+1));
RMU = RMU + ERU^2;
if MAU <= ERU, MAU = ERU; end
end
end
RMU = sqrt(RMU / ((LL-1) * (MM-1)));
% Print errors
% fprintf('T=%6.4f, MAU=%10.4e, RMU=%10.4e\n', (J-1) * K, MAU, RMU);
% Error Tolerance
BIG = 0.0;
for L = 2:LL
for M = 2:MM
FU = abs(U(L, M, J+1) - OLDU(L, M, J+1));
if FU >= BIG, BIG = FU; end
end
end
if BIG <= 10^(-12), break; end
for L = 2:LL
for M = 2:MM
OLDU(L, M, J+1) = U(L, M, J+1);
end
end
end
% Print errors
fprintf('T=%6.4f, KK=%4.0f, MAU=%10.4e, RMU=%10.4e\n', J*K, KK, MAU, RMU);
% Set New Initial Conditions
for L = 1:LL+1
for M = 1:MM+1
U(L, M, J) = U(L, M, J+1);
end
end
end
disp(' ');
toc % CPU time ends
% Plot the results at t=1
for L = 1:LL+1
for M = 1:MM+1
UU(L, M) = U(L, M, JJ+1);
EUU(L, M) = exp(-1) * cos(pi * X(L)) * cos(pi * Y(M));
end
end
figure(1);
colormap summer;
surf(X, Y, UU);
xlabel('X values');
ylabel('Y values');
zlabel('UU values');
title('Numerical Solution at t=1');
figure(2);
colormap winter;
surf(X, Y, EUU);
xlabel('X values');
ylabel('Y values');
zlabel('EUU values');
title('Exact Solution at t=1');