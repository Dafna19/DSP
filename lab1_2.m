% ��� 2�����. u1(t)=sin(2*pi*f1*t), u2(t)=sin(2*pi*f2*t)
clear all %clear workspace
clc %clear command window
close all %�������� ������.����
nfig = 1; %����� ������������ ����
T = 1*(5e-3); % ����� ����.���.������� 1/f1, 1/f2
f1 = 200;
f2 = 1000;
dt = T*1e-3;
t = -T/2 : dt : T/2;
u1 = sin(2*pi*f1*t);
u2 = sin(2*pi*f2*t);
N = length(u1);
% ���������� u1 � u2
figure(nfig); nfig = nfig+1;
plot(t,u1,'-g', 'LineWidth', 2,'DisplayName',['u_1(t)=sin(2\pi' num2str(f1) 't)']);
title('�������');
% ����� ��������� ����: \theta, \pi etc.
hold on
plot(t,u2,'-c', 'LineWidth', 2,'DisplayName',['u_2(t)=sin(2\pi' num2str(f2) 't)']);
grid on
% ��������� ������������
s = scalar(u1,u2,dt,T)
% �����
Nu1 = sqrt(scalar(u1,u1,1,N) )
Nu2 = sqrt(scalar(u2,u2,1,N))
%����� ������ ������ sqrt(1/2) = 0.7

% ��������������� (����)
% ����������������� ����� ?
g1 = u1/Nu1; % ��������������� u1/||u1||
g2 = u2/Nu2;
Ng1 = sqrt(scalar(g1,g1,1,N) ) % ����� g1
Ng2 = sqrt(scalar(g2,g2,1,N))
% g1, g2 - ������������� u1, u2 => ��� ����������������� �����
% => <g1, g2> = 0 (������)
sf = scalar(g1, g2, dt, T)
plot(t,g1,'-y', 'LineWidth', 2,'DisplayName','g_1(t)');
plot(t,g2,'-b', 'LineWidth', 2,'DisplayName','g_2(t)');
legend('show','Orientation','horizontal','Location','NO');
% legend('u_1(t)', 'u_2(t)', 'g_1(t)','g_2(t)') % ����� ���������� ������
% whitebg(gcf,[0.5 0.4 .5]) % ������ ��� RGB; GetCurrentFigure
% 0.2 0.3 0.3; 0.3 0.2 0.3; 0.5 0.4 .5 ������; 0.8 0.5 0.6; 0.57 0.4 .5
% set(gca, 'Color', [0.57 0.5 .5]); %GetCurrentAxes
