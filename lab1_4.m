% ЦОС 4пункт. преобразование Фурье, свойства
clear all %clear workspace
clc %claer command window
close all %закрытие графич.окон

fs = 1000; %частота сигнала
t = 0 : (1/fs)/15 : 4/fs;
u = sin(2*pi*fs*t);
N = length(u);

%матрица спектральных коэф.(для обратного)
F = zeros(N, N);
for k = 0 : N-1
    for n = 0 : N-1
        F(k+1, n+1) = exp(j*2*pi*n*k/N);
    end
end
%матрица спектральных коэф.(для прямого)
FH = conj(F);
F = F./N;   %нормирование
U = u*FH;   %преобразования
u1 = U*F;

plot(t, u, '-m', 'LineWidth', 3, 'DisplayName','u(t)');
title('преобразование Фурье');
hold on
axis([0 2e-3 -1.2 1.2])
grid on
plot(t, u1, ':g', 'LineWidth', 3, 'DisplayName','u1(t)');
hold off

% ЛИНЕЙНОСТЬ

a1 = 5;
a2 = 4;
u2 = sin(2*pi*(fs/3)*t); %второй сигнал
u3 = a1*u + a2*u2;  %суммарный сигнал

U2 = u2*FH; % второй спектр
U4 = a1*U + a2*U2;  %суммарный спектр
u4 = U4*F;  %сигнал для сравнения

figure(2);
plot(t, u, '-b', 'LineWidth', 1, 'DisplayName','u(t)');
title('свойство линейности');
hold on
grid on
plot(t, u2, '-m', 'LineWidth', 1, 'DisplayName','u2(t)');
plot(t, u3, '-r', 'LineWidth', 2, 'DisplayName',['u3(t) = ' num2str(a1) '*u(t)+' num2str(a2) '*u2(t)']);
plot(t, u4, 'pc', 'LineWidth', 2, 'DisplayName','u4(t)');
hold off

% СДВИГ ПО ВРЕМЕНИ

tau = pi; %сдвиг
shifted = t - tau;
u5 = sin(2*pi*fs*shifted); %u = sin(2*pi*fs*t)

U5 = u5*FH; % амплитудный спектр сдвинутого сигнала
FU = atan(imag(U)./real(U));    % фазовый спектр
FU5 = atan(imag(U5)./real(U5));

figure(3);
plot(t, u, '-b', 'LineWidth', 1, 'DisplayName','u(t)');
title('сдвиг по времени');
axis([0 2e-3 -1.2 1.2])
hold on
grid on
plot(t, u5, '-r', 'LineWidth', 1, 'DisplayName','u(t-\tau)');

figure(4);
plot(t, abs(U), '--b', 'LineWidth', 3, 'DisplayName','|U(f)|');
title('амплитудный спектр');
hold on
plot(t, abs(U5), '-r', 'LineWidth', 2, 'DisplayName','|U5(f)| ~ u(t-\tau)');
grid on

figure(5);
%plot(t, FU, '-b', 'LineWidth', 2, 'DisplayName','U(f)');
stem(t, FU, 'DisplayName','U(f)');
title('фазовый спектр');
hold on
% plot(t, FU5, '-r', 'LineWidth', 2, 'DisplayName','U5(f)');
stem(t, FU5, 'r', 'DisplayName','U5(f)');
grid on

% РАВЕНСТВО ПАРСИВАЛЯ

Eu = sum(u.^2)  %сумма значений сигнала в квадрате
EU = sum(abs(U).^2)/N %сумма значений спектра
%{
figure(6); % выводили для проверки
subplot(2,1,1);
plot(t, u, '-m', 'LineWidth', 2, 'DisplayName','u(t)');
subplot(2,1,2);
plot(t, abs(U), '-g', 'LineWidth', 2, 'DisplayName','|U(f)|');
%}


