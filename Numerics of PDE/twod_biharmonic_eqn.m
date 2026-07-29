%----------------------------------------------------- 
% NAME:- KSHITIJ BAHUKHANDI
% South Asian University 
% ENROLLMENT NUMBER:- SAU\AM(M)\2023\06
% 3rd Sem, Numerics of PDE 
% Lab File 09
% Numerical Solution of 2D Biharmonic Equations.
%Uxxxx + 2Uxxyy +Uyyyy = f(x,y) % U(x,y) = sin(pi*x)*sin(pi*y)
% South Asian University.
%--------------------------------------------------------------------------
clear all;
LL= 15;
MM=LL;
H=1.0/(LL+1)
H=1.0/(MM+1);
disp('Numercial solution of 2D biharmonic equation');
disp('         ');
fprintf('LL=%4.2f \n',LL);
disp('         ');
% formation of Nodal points
for L=1:LL+2
    for M=1:MM+2
        X(L)=(L-1)*H;
        Y(M)= (M-1)*H;
    end
end
%exact solution and RHS function
for L=1:LL+2
    for M=1:MM+2
        EU(L,M)=sin(pi*X(L))*sin(pi*Y(M));
        EV(L,M)=-2*pi*pi*EU(L,M);
        F(L,M)=4*pi*pi*pi*pi*EU(L,M);
    end
end
%boundary conditions
for L=1:LL+2
    for M=1:MM+2
        U(L,1)=EU(L,1);
        U(L,MM+2)=EU(L,MM+2);
        U(1,M)=EU(1,M);
        U(LL+2,M)=EU(LL+2,M);
        V(L,1)=EV(L,1);
        V(L,MM+2)=EV(L,MM+2);
        V(1,M)=EV(1,M);
        V(LL+2,M)=EV(LL+2,M);
    end
end
%intial and first approximation
for L=2:LL+1
    for M=2:MM+1
        OLDU(L,M)=0.0;
        OLDV(L,M)=0.0;
        U(L,M)=OLDU(L,M);
        V(L,M)=OLDV(L,M);
    end
end
tic %CPU time starts
%iteration starts
for KK=1:6000;
    for M=2:MM+1
        for L=2:LL+1
            F0=F(L,M);
            F1=F(L+1,M);
            F2=F(L-1,M);
            F3=F(L,M+1);
            F4=F(L,M-1);
            FF=0.5*H*H*(F1+F2+F3+F4+8*F0);
            U0=U(L,M);
            U1=U(L+1,M);
            U2=U(L-1,M);
            U3=U(L,M+1);
            U4=U(L,M-1);
            U5=U(L+1,M+1);
            U6=U(L+1,M-1);
            U7=U(L-1,M+1);
            U8=U(L-1,M-1);
            V0=V(L,M);
            V1=V(L+1,M);
            V2=V(L-1,M);
            V3=V(L,M+1);
            V4=V(L,M-1);
            V5=V(L+1,M+1);
            V6=V(L+1,M-1);
            V7=V(L-1,M+1);
            V8=V(L-1,M-1);
            VV=0.5*H*H*(V1+V2+V3+V4+8*V0);
            U(L,M)=(U5+U6+U7+U8+4*(U1+U2+U3+U4)-VV)/20;
            V(L,M)=(V5+V6+V7+V8+4*(V1+V2+V3+V4)-FF)/20;
        end
    end
    %calculate errosrs
    MAU=0.0;
    MAV=0.0;
    RMU=0.0;
    RMV=0.0;
    for L=2:LL+1
        for M=2:MM+1
            ERU=abs(U(L,M)-EU(L,M));
            ERV=abs(V(L,M)-EV(L,M));
            RMU=RMU+ERU^2;
            RMV=RMV+ERV^2;
            if MAU<=ERU, MAU=ERU; end
            if MAV<=ERV, MAV=ERV; end 
        end
    end
    RMU=sqrt(RMU/ (LL^2));
    RMV=sqrt(RMV/ (LL^2));
    %MAU= max(max(abs(U-EU)));
    %MAU= max(max(abs(U-EU)));
    %print errors
    % fprintf('LL=%4.2f,kk=%4.0f, MAU=%10.4e, MAV=%10.4e, RMU=%10.4e,
    % RMV=%10.4e\n',LL,KK,MAU,MAV,RMU,RMV);
   % fprintf('LL=%4.2f,kk=%4.0f, MAU=%10.4e, MAV=%10.4e, RMU=%10.4e,
    % RMV=%10.4e\n',...,KK,MAU,MAV,RMU,RMV); 
    %Error Tolerance
    %--------------------
    ERR=0.0;
    for L=2:LL+1
        for M=2:MM+1
            FU=abs(U(L,M)-OLDU(L,M));
             if FU>=ERR, ERR=FU; end
        end
    end
    if ERR<= 10^(-12),break; end
    OLDU=U;
    OLDV=V;
end
disp('   ');
    toc % CPU time ends
    figure(1);
    surf(X,Y,U);
    xlabel('X values');
    ylabel('Y values');
    zlabel('U values');
    title('numerical solution of U');
    figure(2);
    surf(X,Y,EU);
    title('exact solution');
    xlabel('X values');
    ylabel('Y values');
    zlabel('EU values');
    figure(3);
    surf(X,Y,V);
    xlabel('X values');
    ylabel('Y values');
    zlabel('V values');
    title('Numerical Solution of V');
    figure(4);
    surf(X,Y,EV);
    title('Exact Solution');
    xlabel('X values');
    ylabel('Y values');
    zlabel('EV values');