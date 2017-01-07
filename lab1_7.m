% ��� 7�����. ����������������� �����������
clear all %clear workspace
clc %clear command window
close all %�������� ������.����
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

%������������� �����������
tmpImg = zeros(rows, ncols);
%����� �����������
newImg = zeros(ncols, nrows); % ����� ���� ��������������
whos newImg

nn = 1 : 1 : cols;
for i = 1:rows
    for t = 1:ncols %������������ �� ������
        tmpImg(i,t) = sum( img(i,:).* sinc(t/n - nn) );
    end
end
% imshow(tmpImg)
whos tmpImg
% �������������� ��� ����.������ sinc
tmpImg = tmpImg';
whos tmpImg

for i = 1:ncols
    for t = 1:nrows %������������ �� ������
        newImg(i,t) = sum( tmpImg(i,:).* sinc(t/m - nn) );
    end
end
whos newImg
% ���������� � ���������� ���������
newImg = newImg';
whos newImg

outIMG = uint8(newImg);

figure(nfig);nfig = nfig+1;
imshow(outIMG);
imwrite(outIMG, 'testmyImg.bmp')



