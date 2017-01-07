% ��� 2.1. �������� �������� SNR
clear all %clear workspace
clc %clear command window
close all %�������� ������.����
nfig = 1;
a = 3;
M = 10000;
X = (2*a).*rand(1,M) - a;

figure(nfig);nfig = nfig+1;
hist(X);
title('����������� ������� X');

% ������������ SNR
SNR = zeros(1,7);
sumX = sum(X.^2);
for R = 1:7
    N = 2^R;
    delta = 2*a/N; % ������� �������� � �����
    %������� �������
    t = zeros(1, N+1);
    for i = 1:N+1
        t(i) = -a + delta*(i-1);
    end
    y = zeros(1, M);
    for i = 1:length(X) % �����������
        for j = 1:N
            if X(i) >= t(j) && X(i) < t(j+1)
                y(i) = (t(j)+t(j+1))/2;
            end
        end
    end
    SNR(R) = 10 * log10( sumX/sum( (X-y).^2 )  );
end
figure(nfig);nfig = nfig+1;
plot(1:7, SNR, 'b*', 'DisplayName', '����������� �������������');
hold on
title('SNR(R)');

% ��������. SNR
SNRt = 10 * log10((2.^(1:7)).^2);
plot(1:7, SNRt, 'r-', 'DisplayName', '�������������');
xlabel('R');
ylabel('SNR');

% ������� �� ����������� ������ N(0, a/3)
sigma = a/3;
Y = randn(1, M)*sigma;
sumY = sum(Y.^2);
for R = 1:7
    N = 2^R;
    delta = 2*a/N; % ������� �������� � �����
    %������� �������
    t = zeros(1, N+1);
    t(1) = -inf;
    t(N+1) = inf;
    for i = 2:N
        t(i) = -a + delta*(i-1);
    end
    y = zeros(1, M);
    for i = 1:length(Y) % �����������
        if Y(i) >= t(1) && Y(i) < t(2)%������� �����
            y(i) = t(2) - delta;
        elseif Y(i) >= t(N) && Y(i) < t(N+1)%������� ������
            y(i) = t(N) + delta;
        else
            for j = 2:N-1
                if Y(i) >= t(j) && Y(i) < t(j+1)
                    y(i) = (t(j)+t(j+1))/2;
                end
            end
        end
    end
    SNR(R) = 10 * log10( sumY/ sum( (Y-y).^2 ) );
end
plot(1:7, SNR, 'bo', 'DisplayName', '������������� �� N(0, \sigma^2)');
legend show

figure(nfig);nfig = nfig+1;
hist(Y);
title('����������� ������� Y');



