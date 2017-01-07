% ЦОС 7пункт. передискретизация изображения
clear all %clear workspace
clc %clear command window
close all %закрытие графич.окон
nfig = 1;
m = 2;
n = 0.5;
IMG = imread('images/airplane.bmp');

figure(nfig);nfig = nfig+1;
imshow(IMG)

img = double(IMG);
rows = size(img,1);
cols = size(img,2);
nrows = m*rows;
ncols = n*cols;

%промежуточное изображение
tmpImg = zeros(rows, ncols);
%новое изображение
newImg = zeros(ncols, nrows); % потом буду переворачивать
whos newImg

nn = 1 : 1 : cols;
for i = 1:rows
    for t = 1:ncols %переделываем по ширине
        tmpImg(i,t) = sum( img(i,:).* sinc(t/n - nn) );
    end
end
% imshow(tmpImg)
whos tmpImg
% переворачиваем для норм.работы sinc
tmpImg = tmpImg';
whos tmpImg

for i = 1:ncols
    for t = 1:nrows %переделываем по высоте
        newImg(i,t) = sum( tmpImg(i,:).* sinc(t/m - nn) );
    end
end
whos newImg
% возвращаем в нормальное состояние
newImg = newImg';
whos newImg

outIMG = uint8(newImg);

figure(nfig);nfig = nfig+1;
imshow(outIMG);
imwrite(outIMG, 'testmyImg.bmp')



