% A-->B-->c
% B-->D

k1 =2;
k2=0.5;
k3=0.3;
h=0.1;
t_final=10;
t_initial=0;

%%EXPLICITING1 PART1
N_explicit= round((t_final-t_initial)/h);
A_explicit =zeros(N_explicit+1,1);
B_explicit =zeros(N_explicit+1,1);
C_explicit =zeros(N_explicit+1,1);
D_explicit =zeros(N_explicit+1,1);
t =zeros(N_explicit+1,1);

%%intial conditions 
A_explicit(1)=1;
B_explicit(1)=0;
C_explicit(1)=0;
D_explicit(1)=0;

for i=1:N_explicit
    t(i+1)=t(i)+h;
    A_explicit(i+1)=A_explicit(i) + h*(-k1*A_explicit(i));
    B_explicit(i+1)=B_explicit(i) + h*(k1*A_explicit(i)-k2*B_explicit(i)-k3*B_explicit(i));
    C_explicit(i+1)=C_explicit(i) + h*(B_explicit(i)*k2);
    D_explicit(i+1)=D_explicit(i) + h*(B_explicit(i)*k3);
end
figure(1)            %graph for explicit (part1)
plot(t,A_explicit,"-r",t,B_explicit,"-b",t,C_explicit,"-g",t,D_explicit,"-m")


%%implicit(real implicit using jacobian )  PART2
N_implicit= ceil((t_final-t_initial)/h);

t=zeros(N_implicit+1,1);

Y=zeros(4,N_implicit+1);
%%intial conditions 
Y(1,1)=1;         %Y(1,:) is showing conc of A at diff i
Y(2,1)=0;         %Y(2,:) is showing conc of b at diff i
Y(3,1)=0;         %Y(3,:) is showing conc of c at diff i   Y----> implicit
Y(4,1)=0;         %Y(4,:) is showing conc of d at diff i 
X=zeros(4,1);
t(1)=0;

for i =2:N_implicit+1
     Y(:,i)=Y(:,i-1);
     t(i)=t(i-1)+h;
    for k=1:10
        j=give_jacobian(h,k1,k2,k3);
        r=give_residual(Y,i,h,k1,k2,k3);
        Y(:,i) = Y(:,i)- inv(j)*(r);
    end 
     
end



figure(2)
plot(t',Y(1,:),"-r")      %graph for implict part2
hold on 
plot(t',Y(2,:),"-b",t',Y(3,:),"-g",t',Y(4,:),"-m")





% %FUNCTIONS
% function R = give_residual(Y,i,h,k1,k2,k3)
% 
%   R=[Y(1,i)-Y(1,i-1)+h*(k1)*Y(1,i);
%      Y(2,i)-Y(2,i-1)-h*(k1*Y(1,i)-k2*Y(2,i)-k3*Y(2,i));  %just show what are the functions of jacobian and residual they are commented here but written in code below to avoid syntax error
%      Y(3,i)-Y(3,i-1)-h*(k2)*Y(2,i);
%      Y(4,i)-Y(4,i-1)-h*k3*Y(2,i)];
% 
% end
% 
% function j = give_jacobian(h,k1,k2,k3)
% 
%  j=[1+(k1)*h,0,0,0;
%     -h*k1,1+k2*h+k3*h,0,0;
%     0,-h*k2,1,0;
%     0,-h*k3,0,1];
% end

%%rk4 part 3
N_RK4=ceil(t_final-t_initial)/h;


y=zeros(4,N_RK4+1);   
y(1,1)=1;           %y(1,:) showing conc of a at different i
y(2,1)=0;           %y(2,:) showing conc of b at different i
y(3,1)=0;           %y(3,:) showing conc of c at different i   y---->rk4
y(4,1)=0;           %y(4,:) showing conc of d at different i

for i=1:N_RK4
    y(1,i+1)=y(1,i)+ h*f_A(i,h,y);
    y(2,i+1)=y(2,i)+ h*f_B(i,h,y);
    y(3,i+1)=y(3,i)+ h*f_C(i,h,y);
    y(4,i+1)=y(4,i)+ h*f_D(i,h,y);
end

figure(3)
plot(t',y(1,:),"-r")     %graph for rk4
hold on 
plot(t',y(2,:),"-b",t',y(3,:),"-g",t',y(4,:),"-m")


%function of A 
function k_A = f_A(i,h,y)
    v1_A= g_A(y(1,i));
    v2_A= g_A(y(1,i)+h/2*v1_A);
    v3_A= g_A(y(1,i)+h/2*v2_A);
    v4_A= g_A(y(1,i)+h/2*v3_A);
    k_A = (v1_A+2*v2_A+2*v3_A+v4_A)/6;
end
function x_A=  g_A(A)
    k1=2;
    x_A = -k1*(A);
end
% function of b
function k_B = f_B(i,h,y)

    v1_B= g_B(y(1,i),y(2,i));
    v2_B= g_B((y(1,i)),(y(2,i)+h/2*v1_B));
    v3_B= g_B((y(1,i)),(y(2,i)+h/2*v2_B));
    v4_B= g_B((y(1,i)),(y(2,i)+h/2*v3_B));
    k_B = (v1_B+2*v2_B+2*v3_B+v4_B)/6;

end
function x_B=  g_B(A,B)
    k1=2;
    k2=0.5;
    k3=0.3;
    x_B = k1*(A)-k2*(B)-k3*(B);
end
%function of c
    function k_C = f_C(i,h,y)
    v1_C= g_C(y(2,i));
    v2_C= g_C(y(2,i));
    v3_C= g_C(y(2,i));
    v4_C= g_C(y(2,i));
    k_C = (v1_C+2*v2_C+2*v3_C+v4_C)/6;
    end


function x_C=  g_C(A)
    k2=0.5;
    x_C= k2*(A);
end
%function of d
function k_D = f_D(i,h,y)
    v1_D= g_D(y(2,i));
    v2_D= g_D(y(2,i));
    v3_D= g_D(y(2,i));
    v4_D= g_D(y(2,i));
    k_D = (v1_D+2*v2_D+2*v3_D+v4_D)/6;
end
function x_D=  g_D(A)
    k3=0.3;
    x_D= k3*(A);
end



%functions for implicit part2
%FUNCTIONS
function R = give_residual(Y,i,h,k1,k2,k3)

  R=[Y(1,i)-Y(1,i-1)+h*(k1)*Y(1,i);
     Y(2,i)-Y(2,i-1)-h*(k1*Y(1,i)-k2*Y(2,i)-k3*Y(2,i));
     Y(3,i)-Y(3,i-1)-h*(k2)*Y(2,i);
     Y(4,i)-Y(4,i-1)-h*k3*Y(2,i)];

end

function j = give_jacobian(h,k1,k2,k3)

 j=[1+(k1)*h,0,0,0;
    -h*k1,1+k2*h+k3*h,0,0;
    0,-h*k2,1,0;
    0,-h*k3,0,1];
end