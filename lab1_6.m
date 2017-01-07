% ЦОС 6пункт. изображения
clear all %clear workspace
clc %claer command window
close all %закрытие графич.окон
nfig = 1;

img1 = 'images/man_first.bmp';
img2 = 'images/man_second.bmp';
A = imread(img1);
B = imread(img2);

AA = A (:,:,1);
BB = B (:,:,1);

for i = 1:length(AA(:,1)) % по строкам
    F1(i,:) = fft(AA(i,:));
    F2(i,:) = fft(BB(i,:));
end
for i = 1:length(F1(1,:)) % по столбцам
    F1(:,i) = fft(F1(:,i));
    F2(:,i) = fft(F2(:,i));
end
F3 = F1.*conj(F2);

for i = 1:length(F3(:,1)) % по столбцам
    C(i,:) = ifft(F3(i,:));
end
for i = 1:length(C(1,:))
    C(:,i) = ifft(C(:,i));
end
figure(nfig);nfig = nfig+1;
surf(real(C))

F3norm = F3./abs(F3);   %нормированное
for i = 1:length(F3norm(:,1)) % обратное
    Cn(i,:) = ifft(F3norm(i,:));
end
for i = 1:length(Cn(1,:))
    Cn(:,i) = ifft(Cn(:,i));
end
Cn = real(Cn);
figure(nfig);nfig = nfig+1;
surf(Cn);
[I,J] = find(Cn == max(max(Cn)))
Cn(I,J) %проверка
corr = zeros(length(Cn(:,1)), length(Cn(1,:)), 3);
corr(I,J,:) = 1;
figure(nfig);nfig = nfig+1;
image(corr);


%{
C = ifftshift(ifftshift(F3),2);% неважно, в каком порядке
figure(nfig); nfig = nfig+1;
image(C)

{
count = 0;
for i = 1:170
    for j = 1:173
        if C(i,j,1) ~= 255
            count = count + 1;
            pxl(1, count) = i;
            pxl(2, count) = j;
        end
    end
end
count
%}



%{
% проверяю БПФ
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
x1 = cos(2*pi*50*t);          % First row wave
x2 = cos(2*pi*150*t);         % Second row wave
x3 = cos(2*pi*300*t);         % Third row wave

X = [x1; x2; x3];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:100),X(i,1:100))
    title(['Row ',num2str(i),' in the Time Domain'])
end

% увеличиваетс размер результата до следующей степени 2
% (было 1000 станет 1024)
n = 2^nextpow2(L);
dim = 2;    % проходимся по строкам
Y = fft(X,n,dim);
P2 = abs(Y/n);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);    % увеличивает амплитуду (с 0,5 до 1)
for i=1:3
    subplot(3,1,i)
    plot(0:(Fs/n):(Fs/2-Fs/n),P1(i,1:n/2))
    title(['Row ',num2str(i), ' in the Frequency Domain'])
end
%}

