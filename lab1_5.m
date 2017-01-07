% ��� 5�����. �����-����
clear all %clear workspace
clc %claer command window
close all %�������� ������.����
nfig = 1;

name = 'DTMF/11.wav';
% FS - ������� �������������
[samples,FS,nbits] = wavread(name);
figure(nfig); nfig = nfig+1;
plot(samples(), 'b');
title(name);
% 10 ������
freq = zeros(10, 2);
for num = 0 : 9
    for i = 1 : 1600
        sam16(i,1) = samples(i + 1600*num,1);
    end
    Samples = abs(fft(sam16));
    Samples = Samples/1600;
    f = linspace(0, FS, 1600);
    
    figure(nfig); nfig = nfig+1;
    plot(f,Samples,'m','DisplayName',name,'LineWidth', 1);
    axis([0, 4000, 0, 0.12])
    title([name '(' num2str(num) ')']);
%}
    [Y,I] = max(Samples); %������� ������ ��������
    freq(num+1,1) = floor(f(I));%������ �������
    if f(I)<1209
        [Y,I1] = max(Samples(I+1:end/2,1) );%������� ������ ��������
        I = I+I1;
    else
        [Y,I] = max(Samples(1:I-1,1) );
    end 
    freq(num+1,2) = floor(f(I));%������ �������

end
freq
%{
1       2       3   	A 	697 ��
4   	5   	6   	B 	770 ��
7       8   	9   	C 	852 ��
*       0   	#   	D 	941 ��
1209   1336    1477   1633 ��
%}




