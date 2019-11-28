function [ output_coher] = coherence( gmaster, gslave )
%UNTITLED1 Summary of this function goes here
%  master slave 分别为主副影像的复数阵
[row,col]=size(gmaster);
output_coher=zeros(row,col);
size_win=3; %相干系统估计窗口大小
U=0;
U1=0;
U2=0;
k = fix(size_win/2);
for i= k+1:row-k
    for j=k+1:col-k
       for n=-k: k
            for m=-k:k
            U= gmaster(i+n,j+m)^2*gslave(i+n,j+m)^2+U;
            U1= gmaster(i+n,j+m)^4+U1;
            U2= gslave(i+n,j+m)^4+U2;
            end
        end
        output_coher(i,j)=U/sqrt(U1*U2);
        if  output_coher(i,j)>0.5
            output_coher(i,j)= sqrt(2*output_coher(i,j)-1);
        else
            output_coher(i,j)=0;
        end
%        if output_coher(i,j) < 0.5    %相干系数阈值选取
 %            output_coher(i,j)=0; 
 %       end
        U=0;
        U1=0;
        U2=0;
    end
end

