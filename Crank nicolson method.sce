//Date : 22/4/2022
//Aim : To solve Heat Equation with Initial Boundary condition using Crank Nicolson method.
clc
clear
clf
//BTCS-backward in time and central in space
q=input("Question to be solved : ")
N=input("Enter No of intervals for space discretizarion= ")
M=input("Enter No of intervals for time discretizarion= ")
t0=0//initial time
tn=1//final time
ht=(tn-t0)/M
disp("del_t = ",ht)
t=t0:ht:tn
select q 
case 1//q1
    a=0;b=1 ;D=1/16
    hx=(b-a)/N
    x=a:hx:b
    fx=2*sin(2*%pi*x)
    ua=zeros(1,M+1)
    ub=zeros(1,M+1)
case 2//q4
    a=0;b=%pi;D=1
    hx=(b-a)/N
    x=a:hx:b
    fx=(cos(x))^2
    ua=(1/2)*(1.+exp(-4*t))
    ub=(1/2)*(1.+exp(-4*t))
case 3//q5
    a=0;b=%pi;D=1
    hx=(b-a)/N
    x=a:hx:b
    fx=(sin(x))^3
    ua=zeros(1,M+1)
    ub=zeros(1,M+1)
case 4//q2
    a=0;b=1;D=1
    hx=(b-a)/N
    x=a:hx:b
    fx=1-x'-(1/%pi)*sin(2*%pi*x')
    ua=ones(1,M+1)
    ub=zeros(1,M+1)
case 5//q3
    a=0;b=%pi;D=1
    hx=(b-a)/N
    x=a:hx:b
    fx=sin(x')
    ua=zeros(1,M+1)
    ub=zeros(1,M+1)
case 6
    a=0;b=1;D=1/25
    hx=(b-a)/N
    x=a:hx:b
    fx=sin(x')
    ua=(-sin(t/5))
    ub=(sin(1-t/5))
end
disp("Coefficient of diffusion D = ",D)
lambda=(1/2)*(D*(ht))/(hx)^2
//disp("The value of lambda is =",lambda)
disp("del_x = ",hx)

//matrix I-
I=eye((N-1),(N-1))

//matrix A-
A=zeros((N-1),(N-1))
for i=1:(N-1)
    A(i,i)=2
    for j=1:(N-1)
        if j==i+1
            A(i,j)=-1
        end
        if j==i-1
            A(i,j)=-1
        end 
    end
end
disp("Tridiagonal matrix = ",A)
E1=(I+lambda*A)
E2=(I-lambda*A)

//omega-
w1=zeros(1,(N-1))
B=zeros(1,(N-1))
B_plus=zeros(1,(N-1))
w1(1,:)=fx(2:N)
omega=zeros((N+1),(M+1))
omega(:,1)=fx'

for i=1:M
    B(1)=-1*ua(i)
    B(N-1)=-1*ub(i)
    B_plus(1,1)=-1*ua(i+1)
    B_plus(1,N-1)=-1*ub(i+1)
    rhs=(E2)*w1'-(lambda*(B'+B_plus'))
    W=inv(E1)*rhs
    omega(:,i+1)=[ua(i+1) W' ub(i+1)]'
    w1=W'
end
for  j=1:M+1
    if q==1 then
        u_exact(:,j)=2*exp(-(((%pi)^2)/4)*t(j))*sin(2*%pi*x')
    elseif q==2
        u_exact(:,j)=(1/2)+(1/2)*(exp(-4*t(j)))*cos(2*x')
    elseif q==3
        u_exact(:,j)=(3/4)*exp(-t(j))*sin(x')-(1/4)*exp(-9*t(j))*sin(3*x')
    elseif q==4
        u_exact(:,j)=1-x'-(1/%pi)*exp(-4*((%pi)^2)*t(j))*sin(2*%pi*x')
    elseif q==5
        u_exact(:,j)=exp(-t(j))*sin(x')
    elseif q==6
        u_exact(:,j)=sin(x')-t(j)/5
    end
end

disp("Solution of heat equation = ",clean(omega))
disp("Exact solution = ",clean(u_exact))

subplot(1,2,1)
plot3d(x,t,omega)
title('Solution by Crank Nikolson method','color','green','Fontsize','5')
subplot(1,2,2)
plot3d(x,t,u_exact)
title('Exact Solution','color','green','Fontsize','5')
