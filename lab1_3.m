% ��� 3�����. ������������� u(t)
clear all %clear workspace
clc %claer command window
close all %�������� ������.����
F = 600;
T = 10/F; %������������ �������
fd1 = 1.5*F; % ������� �������������
fd2 = 1.75*F;
fd3 = 2*F;
fd4 = 3*F;
fd5 = 1000*F;

t = 0 : (1/F)/1000 : T;
u = sin(2*pi*F*t);
plot(t,u,'-m', 'LineWidth', 5,'DisplayName','u(t)');
hold on
axis([-0.0005 0.017 -1.5 1.5])
grid on

Td1 = 1/fd1; % ��� �������������
Td2 = 1/fd2;
Td3 = 1/fd3;
Td4 = 1/fd4;
Td5 = 1/fd5;

% ������� ��������
t1 = 0:Td1:T;
t2 = 0:Td2:T;
t3 = 0:Td3:T;
t4 = 0:Td4:T;
t5 = 0:Td5:T;

u1 = sin(2*pi*F*t1);
u2 = sin(2*pi*F*t2);
u3 = sin(2*pi*F*t3);
u4 = sin(2*pi*F*t4);
u5 = sin(2*pi*F*t5);
% stairs(t1,u1, '-g', 'LineWidth', 2,'DisplayName','u1[n]')

% ��������������
rec1 = zeros(1, length(t));
rec2 = zeros(1, length(t));
rec3 = zeros(1, length(t));
rec4 = zeros(1, length(t));
rec5 = zeros(1, length(t));

for i = 1:length(t)
    w = pi*(t(i) - t1)/Td1; % 1
    buf = sin(w)./w;
    rec1(1,i) = sum(u1.*buf);
    
    w = pi*(t(i) - t2)/Td2; % 2
    buf = sin(w)./w;
    rec2(1,i) = sum(u2.*buf);
    
    w = pi*(t(i) - t3)/Td3; % 3
    buf = sin(w)./w;
    rec3(1,i) = sum(u3.*buf);
    
    w = pi*(t(i) - t4)/Td4; % 4
    buf = sin(w)./w;
    rec4(1,i) = sum(u4.*buf);
    
    w = pi*(t(i) - t5)/Td5; % 5
    buf = sin(w)./w;
    rec5(1,i) = sum(u5.*buf);
end
plot(t,rec1,'-c', 'LineWidth', 2,'DisplayName','u1(t)');
plot(t,rec2,'-r', 'LineWidth', 1,'DisplayName','u2(t)');
plot(t,rec3,'-k', 'LineWidth', 2,'DisplayName','u3(t)');
%� 3 �� ��������� - ��� ������� ����� � ���� � ������
plot(t,rec4,'-y', 'LineWidth', 2,'DisplayName','u4(t)');
plot(t,rec5,'-g', 'LineWidth', 2,'DisplayName','u5(t)');

