%Name= KSHITIJ BAHUKHANDI
%Roll NO = SAU/AM(M)/2023/06
%%% Lab file = 17
%NUMERICAL SOLUTION OF 3D NONLINEAR PARABOLIC EQN WITH A FORCING FUNCTION
%Uxx+Uyy+Uzz=Ut+alpha*U*U+f(x,y,t)
%U(x,y,t)=exp(-t)*cos(pi*x)*cos(pi*y)*cos(pi*z)
%----------------------------------------------------------------------------
clear all;
clc;
lam=0.8;
alpha=10.0;
ll=8;
mm=ll;
nn=ll;
jj=80;
h=1.0/ll;
k=lam*h*h;
disp('NUMERICAL SOLUTION OF 3D NONLINEAR EQUATION');
disp('');
fprintf('LL=%4.2f,K=%6.4f,LAM=%6.2f\n',ll,k,lam);
disp('');

% DIMENSIONS
X=zeros(1,ll+1);
Y=zeros(1,mm+1);
Z=zeros(1,nn+1);
T=zeros(1,jj+1);
TT=zeros(1,jj+1);
U=zeros(ll+1,mm+1,nn+1,jj+1);
EU=zeros(ll+1,mm+1,nn+1,jj+1);
F=zeros(ll+1,mm+1,nn+1,jj+1);
OLDU = zeros(ll+1, mm+1, nn+1, jj+1); % Initialize OLDU

% FORMATION OF NODAL POINTS
for l=1:ll+1
    for m=1:mm+1
        for n=1:nn+1
            for j=1:jj+1
                X(l)=(l-1)*h;
                Y(m)=(m-1)*h;
                Z(n)=(n-1)*h;
                T(j)=(j-1)*k;
                TT(j)=T(j)+0.5*k;
            end
        end
    end
end

% EXACT SOLUTION & RHS FUNCTION
for l=1:ll+1
    for m=1:mm+1
        for n=1:nn+1
            for j=1:jj+1
                EU(l,m,n,j)=exp(-T(j))*cos(pi*X(l))*cos(pi*Y(m))*cos(pi*Z(n));
                EEU(l,m,n,j)=exp(-TT(j))*cos(pi*X(l))*cos(pi*Y(m))*cos(pi*Z(n));
                F(l,m,n,j)=(1-3*pi*pi)*EEU(l,m,n,j)-alpha*EEU(l,m,n,j)^2;
            end
        end
    end
end

% INITIAL CONDITION
for l=1:ll+1
    for m=1:mm+1
        for n=1:nn+1
            U(l,m,n,1)=EU(l,m,n,1);
        end
    end
end

% BOUNDARY CONDITION
for j=1:jj
    for l=1:ll+1
        for m=1:mm+1
            for n=1:nn+1
                U(1,m,n,j+1)=EU(1,m,n,j+1);
                U(ll+1,m,n,j+1)=EU(ll+1,m,n,j+1);
                U(l,1,n,j+1)=EU(l,1,n,j+1);
                U(l,mm+1,n,j+1)=EU(l,mm+1,n,j+1);
                U(l,m,1,j+1)=EU(l,m,1,j+1);
                U(l,m,nn+1,j+1)=EU(l,m,nn+1,j+1);
            end
        end
    end
end

% SOLUTION AT ADVANCED TIME LEVEL
tic % CPU time starts
for j=1:jj
    % INITIAL AND FIRST APPROXIMATION
    for l=2:ll
        for m=2:mm
            for n=2:nn
                OLDU(l,m,n,j+1)=EU(l,m,n,j+1);
            end
        end
    end
    
    for kk=1:6000  % ITERATION STARTS
        % SOLUTION AT INTERNAL GRID POINTS AT EACH TIME LEVEL
        for l=2:ll
            for m=2:mm
                for n=2:nn
                    F0=F(l,m,n,j);
                    R0=U(l,m,n,j);
                    R1=U(l+1,m,n,j);
                    R2=U(l-1,m,n,j);
                    R3=U(l,m+1,n,j);
                    R4=U(l,m-1,n,j);
                    R5=U(l,m,n+1,j);
                    R6=U(l,m,n-1,j);
                    R7=R1+R2-2*R0;
                    R8=R3+R4-2*R0;
                    R9=R5+R6-2*R0;
                    U0=U(l,m,n,j+1);
                    U1=U(l+1,m,n,j+1);
                    U2=U(l-1,m,n,j+1);
                    U3=U(l,m+1,n,j+1);
                    U4=U(l,m-1,n,j+1);
                    U5=U(l,m,n+1,j+1);
                    U6=U(l,m,n-1,j+1);
                    U7=U1+U2;
                    U8=U3+U4;
                    U9=U5+U6;
                    P0=0.25*alpha*k;
                    P1=1+3*lam+2*P0*R0;
                    P2=k*F0+P0*R0*R0-R0-0.5*lam*(R7+R8+R9+U7+U8+U9);
                    AA=P0*U0*U0+P1*U0+P2;
                    DA=2*P0*U0+P1;
                    U(l,m,n,j+1)=U(l,m,n,j+1)-AA/DA;
                end
            end
        end
        
        % CALCULATE ERRORS
        MAU=0.0;
        RMU=0.0;
        for l=2:ll
            for m=2:mm
                for n=2:nn
                    ERU=abs(U(l,m,n,j+1)-EU(l,m,n,j+1));
                    RMU=RMU+ERU^2;
                    if MAU<=ERU ,MAU=ERU;
                    end
                end
            end
        end
        RMU=sqrt(RMU/((ll-1)*(mm-1)*(nn-1)));
        
        % ERROR TOLERANCE
        BIG=0.0;
        for l=2:ll
            for m=2:mm
                for n=2:nn
                    FU=abs(U(l,m,n,j+1)-OLDU(l,m,n,j+1));
                    if FU>=BIG, BIG=FU;end
                end
            end
        end
        if BIG<=10^(-12), break; end
        
        for l=2:ll
            for m=2:mm
                for n=2:nn
                    OLDU(l,m,n,j+1)=U(l,m,n,j+1);
                end
            end
        end
    end
    
    % SET NEW INITIAL CONDITION
    for l=1:ll+1
        for m=1:mm+1
            for n=1:nn+1
                U(l,m,n,j)=U(l,m,n,j+1);
            end
        end
    end
end

% PLOT THE RESULTS AT T=1 AND z=1
for l=1:ll+1
    for m=1:mm+1
        UU(l,m)=U(l,m,nn+1,jj+1);
        EUU(l,m)=exp(-1)*cos(pi*X(l))*cos(pi*Y(m))*cos(pi*Z(nn+1));
    end
end
figure(1);
colormap summer;
surf(X,Y,UU)
xlabel('X values');
ylabel('Y values');
zlabel('UU values');
title('Numerical Solution at t=1')
figure(2);
colormap winter;
surf(X,Y,EUU);
xlabel('X values');
ylabel('Y values');
zlabel('EUU values');
title('Exact Solution at t=1');
