%****************************************************
function [ flag ]=plot_recmap(csar)
%
flag = 1;
N    = 6;
[row,col]=size(csar);
%Ncsar =  zeros(row/N, col);
%for m = 1:row/N
%         Ncsar (m,:) = sum(csar ((m-1)*N+1:1:m*N,:))/N;     
%end
Max=max(max(abs(csar)));
Min=min(min(abs(csar)));
Gsar=(abs(csar)-Min)*256/(Max-Min);
%colormap(gray);
%figure(1);
%imagesc(row,col,Gsar);
c=histeq(Gsar);                  %直方图均衡化
figure;
imshow(c);                   %显示处理后的SAR振幅图像
axis on;
