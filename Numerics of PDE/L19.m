%Name= Kshitij Bahukhandi
%Roll NO.= SAU/AM(M)/2022/06
%LAB:- 19
%SOUTH ASIAN UNIVERSITY
%NUMERICAL SOLUTION OF 3D HEAT CONDUCTION EQUATION with forcing function
%ADI METHOD
%Uxx+Uyy+Uzz=Ut+f(x,y,z,t), U(x,y,t)=sinh(t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
%--------------------------------------------------------------------------
clear all;
clc;
lam=1.6;
ll=8;
mm=ll;
nn=ll;
jj=40;
h=1.0/ll;
k=lam*h*h;
disp('NUMERICAL SOLUTION OF 3D HEAT CONDUCTION EQUATION');
disp('');
fprintf('LL=%4.2f,K=%6.4f,LAM=%6.2f\n',ll,k,lam);
disp('');
%DIMENSIONS
X=zeros(1,ll+1);
Y=zeros(1,mm+1);
Z=zeros(1,nn+1);
T=zeros(1,jj+1);
TT=zeros(1,jj+1);
SUB=zeros(1,ll+1);DIAG=zeros(1,ll+1);
SUP=zeros(1,ll+1);RH=zeros(1,ll+1);
U=zeros(ll+1,mm+1,nn+1,jj+1);EU=zeros(ll+1,mm+1,nn+1,jj+1);
US=zeros(ll+1,mm+1,nn+1);USS=U(ll+1,mm+1,nn+1);
F=zeros(ll+1,mm+1,nn+1,jj+1);EUT=zeros(ll+1,mm+1,nn+1,jj+1);
EEU=zeros(ll+1,mm+1,nn+1,jj+1); EEUT=zeros(ll+1,mm+1,nn+1,jj+1);
%FORMATION OF NODAL POINTS
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
%EXACT SOLUTION & RHS FUNCTION
for l=1:ll+1
    for m=1:mm+1
        for n=1:nn+1
          for j=1:jj+1
            EU(l,m,n,j)=sinh(T(j))*sin(pi*X(l))*sin(pi*Y(m))*sin(pi*Z(n));
            EUT(l,m,n,j)=cosh(T(j))*sin(pi*X(l))*sin(pi*Y(m))*sin(pi*Z(n));
            EEU(l,m,n,j)=sinh(TT(j))*sin(pi*X(l))*sin(pi*Y(m))*sin(pi*Z(n));
            EEUT(l,m,n,j)=cosh(TT(j))*sin(pi*X(l))*sin(pi*Y(m)*sin(pi*Z(n)));
            F(l,m,n,j)=-3*pi*pi*EEU(l,m,n,j)-EEUT(l,m,n,j);
          end
        end
    end
end
%INTIAL CONDITION
for l=1:ll+1
    for m=1:mm+1
        for n=1:nn+1
           U(l,m,n,1)=EU(l,m,n,1);
        end
    end
end
%BOUNDARY CONDITION
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
%SOLUTION AT ADVANCED TIME LEVEL
tic %CPU time starts
for j=1:jj
    %D'YAKONOV ADI  METHOD
    %SECONDARY BOUNDARY CONDTION FOR US
    for l=2:ll
        for m=1:mm:mm+1
            for n=2:nn
                US(l,m,n)=-0.5*lam*U(l-1,m,n,j+1)+(1+lam)*U(l,m,n,j+1)-0.5*lam*U(l+1,m,n,j+1);
            end
        end 
    end
    %first intermediate boundary condtion for uss
    for l=2:ll
        for m=2:mm
            for n=1:nn:nn+1
                USS(l,m,n)=-0.5*lam*US(l,m-1,n)+(1+lam)*US(l,m,n)-0.5*lam*US(l,m+1,n);
            end
        end
    end
    %SOLUTION AT FIRST INTERMEDIATE STEP
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
                R7=U(l+1,m+1,n,j);
                R8=U(l+1,m-1,n,j);
                R9=U(l-1,m+1,n,j);
                R10=U(l-1,m-1,n,j);
                R11=U(l,m+1,n+1,j);
                R12=U(l,m+1,n-1,j);
                R13=U(l,m-1,n+1,j);
                R14=U(l,m-1,n-1,j);
                R15=U(l+1,m,n+1,j);
                R16=U(l+1,m,n-1,j);
                R17=U(l-1,m,n+1,j);
                R18=U(l-1,m,n-1,j);
                R19=U(l+1,m+1,n+1,j);
                R20=U(l+1,m-1,n+1,j);
                R21=U(l-1,m+1,n+1,j);
                R22=U(l-1,m-1,n+1,j);
                R23=U(l+1,m+1,n-1,j);
                R24=U(l+1,m-1,n-1,j);
                R25=U(l-1,m+1,n-1,j);
                R26=U(l-1,m-1,n-1,j);
                R31=R1+R2;
                R32=R3+R4;
                R33=R5+R6;
                R34=R31-2*R0;
                R35=R32-2*R0;
                R36=R33-2*R0;
                R37=R7+R8+R9+R10;
                R38=R11+R12+R13+R14;
                R39=R15+R16+R17+R18;
                R40=R19+R20+R21+R22+R23+R24+R25+R26;
                R41=R37-2*(R31+R32)+4*R0;
                R42=R38-2*(R32+R33)+4*R0;
                R43=R39-2*(R31+R33)+4*R0;
                R44=R40-2*(R37+R38+R39)+4*(R31+R32+R33)-8*R0;
                SUB(n)=-0.5*lam;
                DIAG(n)=1+lam;
                SUP(n)=-0.5*lam;
                RH(n)=R0+0.5*lam*(R34+R35+R36)+0.25*lam*lam*(R41+R42+R43)+0.125*lam*lam*lam*R44-k*F0;
            end
            RH(2)=RH(2)-SUB(2)*USS(l,m,1);
            RH(nn)=RH(nn)-SUP(nn)*USS(l,m,nn+1);
            SUB(2)=0.0;
            SUP(nn)=0.0;
            %TRI DIAGONAL SOLVER
            SUP(2)=SUP(2)/DIAG(2);
            RH(2)=RH(2)/DIAG(2);
            for n=3:nn-1
                SUP(n)=SUP(n)/(DIAG(n)-SUB(n)*SUP(n-1));
            end
            for n=3:nn
                RH(n)=(RH(n)-SUB(n)*RH(n-1))/(DIAG(n)-SUB(n)*SUP(n-1));
            end
            USS(l,m,nn)=RH(nn);
            for n=2:nn-1
                USS(l,m,nn-n+1)=RH(nn-n+1)-SUP(nn-n+1)*USS(l,m,nn-n+2);
            end
        end
    end
    %SOLUTION AT SECOND INTERMEDIATE STEP
    for l=2:ll
        for n=2:nn
            for m=2:mm
                SUB(m)=-0.5*lam;
                DIAG(m)=1+lam;
                SUP(m)=-0.5*lam;
                RH(m)=USS(l,m,n);
            end
            RH(2)=RH(2)-SUB(2)*US(l,1,n);
            RH(mm)=RH(mm)-SUP(mm)*US(l,mm+1,n);
            SUB(2)=0.0;
            SUP(mm)=0.0;
            %TRI DIAGONAL SOLVER
            SUP(2)=SUP(2)/DIAG(2);
            RH(2)=RH(2)/DIAG(2);
            for m=3:mm-1
                SUP(m)=SUP(m)/(DIAG(m)-SUB(m)*SUP(m-1));
            end
            for m=3:mm
                RH(m)=(RH(m)-SUB(m)*RH(m-1))/(DIAG(m)-SUB(m)*SUP(m-1));
            end
            US(l,mm,n)=RH(mm);
            for m=2:mm-1
                US(l,mm-m+1,n)=RH(mm-m+1)-SUP(mm-m+1)*US(l,mm-m+2,n);
            end
        end
    end
    %SOLUTION AT MAIN TIME STEP
    for m=2:mm
        for n=2:nn
            for l=2:ll
                SUB(l)=-0.5*lam;
                DIAG(l)=1+lam;
                SUP(l)=-0.5*lam;
                RH(l)=US(l,m,n);
            end
            RH(2)=RH(2)-SUB(2)*U(1,m,n,j+1);
            RH(ll)=RH(ll)-SUP(ll)*U(ll+1,m,n,j+1);
            SUB(2)=0.0;
            SUP(ll)=0.0;
            %TRI DIAGONAL SOLVER
            SUP(2)=SUP(2)/DIAG(2);
            RH(2)=RH(2)/DIAG(2);
            for l=3:ll-1
                SUP(l)=SUP(l)/(DIAG(l)-SUB(l)*SUP(l-1));
            end
            for l=3:ll
                RH(l)=(RH(l)-SUB(l)*RH(l-1))/(DIAG(l)-SUB(l)*SUP(l-1));
            end
            U(ll,m,n,j+1)=RH(ll);
            for l=2:ll-1
                U(ll-l+1,m,n,j+1)=RH(ll-l+1)-SUP(ll-l+1)*U(ll-l+2,m,n,j+1);
            end
        end
    end
    %CALCULATE ERRORS
    MAU=0.0;
    RMU=0.0;
    for l=2:ll
        for m=2:mm
            for n=2:nn
                ERU=abs(U(l,m,n,j+1)-EU(l,m,n,j+1));
                RMU=RMU+ERU^2;
                if MAU<=ERU, MAU=ERU;
                end
            end
        end
    end
    RMU=sqrt(RMU/((ll-1)*(mm-1)*(nn-1)));
    %PRINT ERRORS
    fprintf('T=%6.4f,MAU=%10.4e,RMU=%10.4e\n',j*k,MAU,RMU);
    %SET NEW INITIAL CONDITION
    for l=1:ll+1
        for m=1:mm+1
            for n=1:nn+1
                U(l,m,n,j)=U(l,m,n,j+1);
            end
        end
    end
end
disp('');
toc %CPU TIME END