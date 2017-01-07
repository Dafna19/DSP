% ЦОС 2.2. анализ скалярного квантования
clear all %clear workspace
clc %clear command window
close all %закрытие графич.окон
nfig = 1;
IMG1 = imread('images/Bright_2.bmp');
IMG2 = imread('images/baboon.ppm');
IMG3 = imread('images/dark_2.bmp');
%перевод в Y
img1 = 0.299*IMG1(:,:,1) + 0.587*IMG1(:,:,2) + 0.114*IMG1(:,:,3);
img = img1';
img11 = double(img(:)');%длинные массивы

img2 = 0.299*IMG2(:,:,1) + 0.587*IMG2(:,:,2) + 0.114*IMG2(:,:,3);
img = img2';
img22 = double(img(:)');

img3 = 0.299*IMG3(:,:,1) + 0.587*IMG3(:,:,2) + 0.114*IMG3(:,:,3);
img = img3';
img33 = double(img(:)');

figure(nfig);nfig = nfig+1;
imhist(img1); 
title('гистограмма засвеченного изображения');
figure(nfig);nfig = nfig+1;
imhist(img2);
title('гистограмма сбалансированного изображения');
hold on
% figure(nfig);nfig = nfig+1;
% imhist(img3);
% title('гистограмма затемненного изображения');

SNR1 = zeros(1,7);  %1
PSNR1 = zeros(1,7);
sum1 = sum(img11.^2);
H1 = zeros(1,7);
SNR2 = zeros(1,7);  %2
PSNR2 = zeros(1,7);
sum2 = sum(img22.^2);
H2 = zeros(1,7);
SNR3 = zeros(1,7);  %3
PSNR3 = zeros(1,7);
sum3 = sum(img33.^2);
H3 = zeros(1,7);

% равноменрное скалярное квантование
%{
for R = 1:7
    N = 2^R;%кол-во квантов
    delta = 255/N;
    t = zeros(1,N+1);
    for i = 1:N+1
        t(i) = delta*(i-1);
    end
    y1 = zeros(1,length(img11)); %проквантованные
    for i = 1: length(img11)
        for j = 1:N
            if img11(i) >= t(j) && img11(i) < t(j+1)
                y1(i) = (t(j)+t(j+1))/2;
            end
        end        
    end
    SNR1(R) = 10 * log10( sum1/sum((img11-y1).^2) );
    PSNR1(R) = 10* log10(350*467*255^2 / sum((img11-y1).^2) );
    %kkk=sum((img11-y1).^2)/length(img11)
    
    y2 = zeros(1,length(img22));
    for i = 1: length(img22)
        for j = 1:N
            if img22(i) >= t(j) && img22(i) < t(j+1)
                y2(i) = (t(j)+t(j+1))/2;
            end
        end        
    end
    SNR2(R) = 10 * log10( sum2/sum((img22-y2).^2) );
    PSNR2(R) = 10* log10(480*500*255^2 / sum((img22-y2).^2) );
    
    y3 = zeros(1,length(img33));
    for i = 1: length(img33)
        for j = 1:N
            if img33(i) >= t(j) && img33(i) < t(j+1)
                y3(i) = (t(j)+t(j+1))/2;
            end
        end        
    end
    SNR3(R) = 10 * log10( sum3/sum((img33-y3).^2) );
    PSNR3(R) = 10* log10(486*648*255^2 / sum((img33-y3).^2) );
    
    %разлёт заквантованных значений 
    x = delta*0.5 : delta : (2^R-1)*delta+0.5*delta;
    p = hist(y1, x); %p - количество значений
    p(~p)=[]; % удаляет нулевые элементы
    summa = sum(p);
    p = p./summa; %делаем вероятности
    H1(R) = -sum( p.*log2(p));
    % для y2
    p = hist(y2, x);
    p(~p)=[];
    summa = sum(p);
    p = p./summa;
    H2(R) = -sum( p.*log2(p));
    % для y3
    p = hist(y3, x);
    p(~p)=[];
    summa = sum(p);
    p = p./summa;
    H3(R) = -sum( p.*log2(p));
end
figure(nfig);nfig = nfig+1;
plot(1:7, SNR1, 'r*-', 'DisplayName', 'светлое');
hold on
title('SNR(R) (равномерное с.к.)');
plot(1:7, SNR2, 'b*-', 'DisplayName', 'нормальное');
plot(1:7, SNR3, 'k*-', 'DisplayName', 'темное');

figure(nfig);nfig = nfig+1;
plot(1:7, PSNR1, 'ro-', 'DisplayName', 'светлое');
hold on
title('PSNR(R) (равномерное с.к.)');
plot(1:7, PSNR2, 'bo-', 'DisplayName', 'нормальное');
plot(1:7, PSNR3, 'ko-', 'DisplayName', 'темное');

figure(nfig);nfig = nfig+1;
plot(H1,PSNR1, 'rs:', 'DisplayName', 'светлое');
hold on
title('PSNR(H) (равномерное с.к.)');
plot(H2, PSNR2, 'bs:', 'DisplayName', 'нормальное');
plot(H3, PSNR3, 'ks:', 'DisplayName', 'темное');
%}
% неравноменрное скалярное квантование
%{
% светлое изображение
errors = [209 70 15 3 0.5 0.2 0.07]; % для Ts
deltas = [1 1 1 0.1 0.01 0.0002 0.01]; %для deltaS 
for R = 1:7
    N = 2^R; % число квантов
    % k = 1; % номер итерации
    %%%%%%%%%%%%%%%%%  ПОДОБРАТЬ?????
    Tk = 10; % max число итераций
    deltaS = deltas(R); % порог для разности ошибок квантования
    Ts = errors(R); % порог для ошибки квантования
    %%%%%%%%%%%%%%%%%%
    prevEx = 1000; % предыдущая ошибка квантования
    t = zeros(1, N+1); %границы для каждой итерации
    
    delta = 255/N;
    for i = 1:N
        t(i) = delta*(i-1); % границы как для равномерного квантования
    end
    t(N+1) = inf;
    %была инициализация
    
    for k = 1:Tk
        [Q,ind] = histc(img11,t);
        Q(end) = []; % убираем последний столбец,
        %который отвечает за числа=inf  и который сам равен 0
        
        %проверяем, что кванты не пусты
        test = 0;
        while test == 0
            test = 1;
            shift = 0; % сдвиг для индекса, с которого начнем менять границы
            for i = length(Q):-1:1
                if Q(i) == 0
                    test = 0;
                    shift = shift + 1;
                    %раздвигаем границы вправо
                    t(end) = 255;
                    for j = length(t)-shift : -1 : i+1
                        t(j) = floor( (t(j) + t(j+1))/2); %floor чтобы не наезжали друг на друга
                    end
                end
            end
            if test == 0 % если нули были
                t(end) = inf;
                [Q,ind] = histc(img11,t);
                Q(end) = [];
            end
        end
        % nulls = zeros(1,length(t)); t(end)=255; plot(t, nulls,'r*','DisplayName','границы квантов')
        
        % шаг 1
        y = zeros(1, length(Q)); % средние значения
        for i = 1:length(img11)
            y(ind(i)) = y(ind(i)) + img11(i);
        end
        for i = 1:length(y)
            y(i) = y(i)/Q(i); % формула (3.20)
        end
        
        % шаг 2
        Ex = sum( (img11 - y(ind)).^2)/length(img11) % (3.21)
        
        %шаг 3
        if Ex < Ts || abs(Ex-prevEx)<deltaS || k == Tk
            break;%ломает k
        end
        prevEx = Ex;
        
        % шаг 4
        for i = 2:length(t)-1
            t(i) = (y(i)+y(i-1))/2; % (3.16)
        end
        % шаг 5 k++
    end    
    PSNR1(R) = 10* log10(350*467*255^2 / sum((img11-y(ind)).^2) ); % (3.26)
    summa = sum(Q);
    p = Q./summa;
    H1(R) = -sum( p.*log2(p));
end
%
% темное изображение
errors = [100 50 7 2.5 0.5 0.11 0.0001]; % для Ts
deltas = [10 1 1 0.5 0.05 0.001 0.00001]; %для deltaS 
for R = 7:7
    N = 2^R; % число квантов
    % k = 1; % номер итерации
    %%%%%%%%%%%%%%%%%  ПОДОБРАТЬ?????
    Tk = 10; % max число итераций
    deltaS = deltas(R); % порог для разности ошибок квантования
    Ts = errors(R); % порог для ошибки квантования
    %%%%%%%%%%%%%%%%%%
    prevEx = 1000; % предыдущая ошибка квантования
    t = zeros(1, N+1); %границы для каждой итерации
    
    delta = 255/N;
    for i = 1:N
        t(i) = delta*(i-1); % границы как для равномерного квантования
    end
    t(N+1) = inf;
    %была инициализация
    
    for k = 1:Tk
        [Q,ind] = histc(img33,t);
        Q(end) = []; % убираем последний столбец,
        %который отвечает за числа=inf  и который сам равен 0
        
        %проверяем, что кванты не пусты
        test = 0;
        counts = 0;
        while test == 0
            test = 1;
            shift = 0; % сдвиг для индекса, с которого начнем менять границы
            for i = 1:length(Q)
                if Q(i) == 0
                    test = 0;
                    shift = shift + 1;
                    %раздвигаем границы влево
                    t(end) = 255;
                    for j = 1+shift : i
                        t(j) = ceil( (t(j) + t(j-1))/2); %ceil чтобы не наезжали друг на друга
                    end
                end
            end
            if test == 0 % если нули были
                counts = counts + 1;
                t(end) = inf;
                [Q,ind] = histc(img33,t);
                Q(end) = [];
            end
            if counts >30
                Q(~Q)=0.1;
                break; %если не получается убрать нулевые кванты(при 128 не получается)
            end
        end
        % nulls = zeros(1,length(t)); t(end)=255; plot(t, nulls,'r*','DisplayName','границы квантов')
        
        % шаг 1
        y = zeros(1, length(Q)); % средние значения
        for i = 1:length(img33)
            y(ind(i)) = y(ind(i)) + img33(i);
        end
        for i = 1:length(y)
            y(i) = y(i)/Q(i); % формула (3.20)
        end
        
        % шаг 2
        Ex = sum( (img33 - y(ind)).^2)/length(img33) % (3.21)
        
        %шаг 3
        if Ex < Ts || abs(Ex-prevEx)<deltaS || k == Tk
            break;%ломает k
        end
        prevEx = Ex;
        
        % шаг 4
        for i = 2:length(t)-1
            t(i) = (y(i)+y(i-1))/2; % (3.16)
        end
        % шаг 5 k++
    end    
    PSNR3(R) = 10* log10(486*648*255^2 / sum((img33-y(ind)).^2) ); % (3.26)
    summa = sum(Q);
    Q(~Q)=[]; % в случае 128
    p = Q./summa;
    H3(R) = -sum( p.*log2(p));
end
%}
% сбалансированное изображение
%
errors = [600 160 40 15 5]; % для Ts
deltas = [10 5 2 2 1]; %для deltaS 
middle = (min(img22) + max(img22))/3; %"середина" значений
for R = 6:6
    N = 2^R; % число квантов
    % k = 1; % номер итерации
    %%%%%%%%%%%%%%%%%  ПОДОБРАТЬ?????
    Tk = 10; % max число итераций
    deltaS = 0.1; % порог для разности ошибок квантования
    Ts = 0.5; % порог для ошибки квантования
    %%%%%%%%%%%%%%%%%%
    prevEx = 1000; % предыдущая ошибка квантования
    t = zeros(1, N+1); %границы для каждой итерации
    
    delta = 255/N;
    for i = 1:N
        t(i) = delta*(i-1); % границы как для равномерного квантования
    end
    t(N+1) = inf;
    %была инициализация
    
    for k = 1:Tk
        [Q,ind] = histc(img22,t);
        Q(end) = []; % убираем последний столбец,
        %который отвечает за числа=inf  и который сам равен 0
        
        %проверяем, что кванты не пусты
        test = 0;
        counts = 0;
        while test == 0
            test = 1;
            shift1 = 5;shift2 = 5; % сдвиг для индекса, с которого начнем менять границы
            for i = length(Q):-1:1
                if Q(i) == 0
                    test = 0;
                    index = find(t > middle,1); %первый индекс t, которое > middle
                    t(end) = 255;
                    if t(i) >= middle
                        shift1 = shift1 + 1;
                        %раздвигаем границы влево
                        for j = index+shift1 : i
                            t(j) = ceil( (t(j) + t(j-1))/2);
                        end
                    else
                        shift2 = shift2 + 1;
                        %раздвигаем границы вправо
                        for j = index-1-shift2 : -1 : i+1
                            t(j) = floor( (t(j) + t(j+1))/2); %floor чтобы не наезжали друг на друга
                        end
                    end
                end
            end
            if test == 0 % если нули были
                counts = counts + 1;
                t(end) = inf;
                [Q,ind] = histc(img22,t);
                Q(end) = [];
            end
            if counts >30
                Q(~Q) = 0.1;
                break; %если не получается убрать нулевые кванты
            end
        end
         nulls = zeros(1,length(t)); t(end)=255; plot(t, nulls,'r*','DisplayName','границы квантов')
        
        % шаг 1
        y = zeros(1, length(Q)); % средние значения
        for i = 1:length(img22)
            y(ind(i)) = y(ind(i)) + img22(i);
        end
        for i = 1:length(y)
            y(i) = y(i)/Q(i); % формула (3.20)
        end
        if counts>30 %экстренный случай
            for i = length(y)-1:-1:2
                if y(i)==0 && (y(i-1)~=0 || y(i+1)~=0 )
                    y(i)=(y(i-1)+y(i+1)+1)/2;
                end
            end
        end
        
        % шаг 2
        Ex = sum( (img22 - y(ind)).^2)/length(img22) % (3.21)
        
        %шаг 3
        if Ex < Ts || abs(Ex-prevEx)<deltaS || k == Tk
            break;%ломает k
        end
        prevEx = Ex;
        
        % шаг 4
        for i = 2:length(t)-1
            t(i) = (y(i)+y(i-1))/2; % (3.16)
        end
        % шаг 5 k++
    end    
    PSNR2(R) = 10* log10(480*500*255^2 / sum((img22-y(ind)).^2) ); % (3.26)
    summa = sum(Q);
    p = Q./summa;
    H2(R) = -sum( p.*log2(p));
end

figure(nfig);nfig = nfig+1;
plot(1:7, PSNR1, 'ro-', 'DisplayName', 'светлое');
title('PSNR(R) (неравномерное с.к.)');
hold on
plot(1:7, PSNR3, 'ko-', 'DisplayName', 'темное');
plot(1:7, PSNR2, 'bo-', 'DisplayName', 'нормальное');

figure(nfig);nfig = nfig+1;
plot(H1,PSNR1, 'rs:', 'DisplayName', 'светлое');
title('PSNR(H) (неравномерное с.к.)');
hold on
plot(H3, PSNR3, 'ks:', 'DisplayName', 'темное');
plot(H2, PSNR2, 'bs:', 'DisplayName', 'нормальное');
%}

% (для img11)
% R Ex k=1    k=2        k=3       k=4       k=5      k=6     k=7     k=8
% 1 791.12    296.6263   218.2732  210.2088  209.3253
% 2 98.87     91.4216    86.2135   82.0144   79.0022  76.1425 73.2244 70.9431
% (2)------ k=9 69.1727
% 3 30.85     22.5343    19.9053   19.06
% 4 4.74      4.5054     4.4487 (0.1)
% 5 1.1215    1.1113     1.1069 (0.01)
% 6 0.2966    0.2963     0.2963
% 7 0.0858    0.0858
% для img33 ----------------------------
% 1 413.5958  369.4467  92.0266
% 2 57.4749   52.3004   51.9860
% 3 8.5175    7.4482    7.1995
% 4 2.815     2.6331
% 5 0.5644    0.527
% 6 0.1167    0.1166
% 7 3e-5
% для img22 -----------------------------
% 1 715.993   644.3     620.81  615.82
% 2 259.58    217.17    202.63  188.08  179.94  172.8   167.96
% 3 63.49     53.56     48.89
% 4 16.89     15.003
% 5 6.065     5.0081    4.82
% 6 2.2581
% 7 


ab=[0 11 0 25 63 0 0 223 255];
for i=length(ab)-1:-1:2
    if ab(i)==0 && (ab(i-1)~=0 || ab(i+1)~=0 )
        ab(i)=(ab(i-1)+ab(i+1)+1)/2;
    end
end
%i=find(ab>120,1)



