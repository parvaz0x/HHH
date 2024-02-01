
clc;
clear all;
close all;
tic;
% Cost function weight matrices, Eq. (2)
R = eye(2);
Q= diag([10^4,10^6]);
F=diag([1,1]);

mg=0.46;
mx=4.9;
my=8.5;
r=0.12;
J=0.05;
d=1.2;
g=9.81;

t0=0;
tf=20;
DeltaT=0.01;
t=0:DeltaT:9.999;

it= round((tf-t0)/DeltaT); %Number of steps during each simulation 

% Memory allocation
x = zeros(6,it); % History of states
u = zeros(2,it); % History of control

x(:,1) =  [0 0 -pi/4 0 0 0]'; % Initial state vector/ Initial condition;
%% Online Control using Algorithm 1
% Euler integration is used for propagating the continuous dynamics, given
% sampling time DeltaT

for i = 2:it 
    % Step 1
    % Updating the state dependent matrices
   Ax=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1
    0 0 -mg*g*sin(x(3,i-1))/(mx*x(3,i-1)) -d/mx 0 0
    0 0 mg*g*(cos(x(3,i-1))-1)/(my*x(3,i-1)) 0 -d/my 0
    0 0 0 0 0 0];

    Bx=[0 0;0 0;0 0;cos(x(3,i-1))/mx -sin(x(3,i-1))/mx;
    sin(x(3,i-1))/my cos(x(3,i-1))/my;r/J 0];
    C=[1 0 0 0 0 0;0 1 0 0 0 0];
    
    t(i)=(i)*DeltaT;
    
    Z=[0.5*t(i);0.1*sin(t(i))]; 
    
    % Step 2
    % Solving Eq. (41) for its negative definite solution
    Pss = -care(-Ax,Bx,C'*Q*C,R);
    Acl = Ax - Bx*R^-1*Bx'*Pss;
    D = lyap(Acl,-Bx*R^-1*Bx');
    K_tf = inv(C'*F*C - Pss);
    g_tf=C'*F*Z;
    K = expm(Acl*(t(i)-tf))*(K_tf - D)*expm(Acl'*(t(i)-tf)) + D;
    P = Pss + K^-1;     
    u(:,i) = (-R^-1)*Bx'*(P*x(:,i-1));
    u=min(u,+50);
    u=max(u,-50);
    
    % Applying the control and propagating the states
    x(:,i) = x(:,i-1) + DeltaT*(Ax*x(:,i-1)+Bx*u(:,i-1));
     
end

fprintf('Done \n')

%% Plotting the results

figure;
plot(t,x(1,1:it),'r','LineWidth',2);
xlabel('time(s)');
ylabel('The First State (x)');
grid on;

figure;
plot(t,x(2,1:it),'b','LineWidth',2);
xlabel('time(s)');
ylabel('The second state (y)');
grid on;

figure;
plot(t,x(3,1:it),'b','LineWidth',2);
xlabel('time(s)');
ylabel('pitch Angle(Teta)');
grid on;

figure;
plot(t,u(1,1:it),'r','LineWidth',2);
hold on;
plot(t,u(2,1:it),'LineWidth',2);
legend('u_1','u_2')
title('Control');
xlabel('time(s)');
ylabel('u');
grid on;
toc
