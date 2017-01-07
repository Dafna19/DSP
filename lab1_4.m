% ��� 4�����. �������������� �����, ��������
clear all %clear workspace
clc %claer command window
close all %�������� ������.����

fs = 1000; %������� �������
t = 0 : (1/fs)/15 : 4/fs;
u = sin(2*pi*fs*t);
N = length(u);

%������� ������������ ����.(��� ���������)
F = zeros(N, N);
for k = 0 : N-1
    for n = 0 : N-1
        F(k+1, n+1) = exp(j*2*pi*n*k/N);
    end
end
%������� ������������ ����.(��� �������)
FH = conj(F);
F = F./N;   %������������
U = u*FH;   %��������������
u1 = U*F;

plot(t, u, '-m', 'LineWidth', 3, 'DisplayName','u(t)');
title('�������������� �����');
hold on
axis([0 2e-3 -1.2 1.2])
grid on
plot(t, u1, ':g', 'LineWidth', 3, 'DisplayName','u1(t)');
hold off

% ����������

a1 = 5;
a2 = 4;
u2 = sin(2*pi*(fs/3)*t); %������ ������
u3 = a1*u + a2*u2;  %��������� ������

U2 = u2*FH; % ������ ������
U4 = a1*U + a2*U2;  %��������� ������
u4 = U4*F;  %������ ��� ���������

figure(2);
plot(t, u, '-b', 'LineWidth', 1, 'DisplayName','u(t)');
title('�������� ����������');
hold on
grid on
plot(t, u2, '-m', 'LineWidth', 1, 'DisplayName','u2(t)');
plot(t, u3, '-r', 'LineWidth', 2, 'DisplayName',['u3(t) = ' num2str(a1) '*u(t)+' num2str(a2) '*u2(t)']);
plot(t, u4, 'pc', 'LineWidth', 2, 'DisplayName','u4(t)');
hold off

% ����� �� �������

tau = pi; %�����
shifted = t - tau;
u5 = sin(2*pi*fs*shifted); %u = sin(2*pi*fs*t)

U5 = u5*FH; % ����������� ������ ���������� �������
FU = atan(imag(U)./real(U));    % ������� ������
FU5 = atan(imag(U5)./real(U5));

figure(3);
plot(t, u, '-b', 'LineWidth', 1, 'DisplayName','u(t)');
title('����� �� �������');
axis([0 2e-3 -1.2 1.2])
hold on
grid on
plot(t, u5, '-r', 'LineWidth', 1, 'DisplayName','u(t-\tau)');

figure(4);
plot(t, abs(U), '--b', 'LineWidth', 3, 'DisplayName','|U(f)|');
title('����������� ������');
hold on
plot(t, abs(U5), '-r', 'LineWidth', 2, 'DisplayName','|U5(f)| ~ u(t-\tau)');
grid on

figure(5);
%plot(t, FU, '-b', 'LineWidth', 2, 'DisplayName','U(f)');
stem(t, FU, 'DisplayName','U(f)');
title('������� ������');
hold on
% plot(t, FU5, '-r', 'LineWidth', 2, 'DisplayName','U5(f)');
stem(t, FU5, 'r', 'DisplayName','U5(f)');
grid on

% ��������� ���������

Eu = sum(u.^2)  %����� �������� ������� � ��������
EU = sum(abs(U).^2)/N %����� �������� �������
%{
figure(6); % �������� ��� ��������
subplot(2,1,1);
plot(t, u, '-m', 'LineWidth', 2, 'DisplayName','u(t)');
subplot(2,1,2);
plot(t, abs(U), '-g', 'LineWidth', 2, 'DisplayName','|U(f)|');
%}


