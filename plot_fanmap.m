function [ flag ] = plot_fanmap( theta,array_map)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
[row,col]=size(array_map);
Min_t=(180-theta)/2;
Max_t=Min_t+theta;
t=linspace(Min_t,Max_t,col); 
Max_r=row/2;
X=zeros(row,col);
Y=zeros(row,col);
Z=zeros(row,col);
for i=1:row
     for  j=1:col
     X(i,j)=i/2*cos(t(j)*pi/180);
      Y(i,j)=i/2*sin(t(j)*pi/180);       
   %    plot(x,y,'y-');
    %   hold on;
     end
end
figure;
C=fliplr(flipud(array_map));
%Max=max(max(abs(C)));
%Min=min(min(abs(C)));
%C=(abs(C)-Min)*256/(Max-Min); %线性拉伸
%C=histeq(C);                  %直方图均衡化

mesh(X,Y,Z,C,'LineWidth',2);
alpha(0.5);
set(gcf,'color','white'); %图形背景设为白色
X_axis=(Max_r)*cos(Min_t*pi/180);
axis([-X_axis,X_axis,0,Max_r]);
grid on; colorbar;
%hold off;
flag=1;
