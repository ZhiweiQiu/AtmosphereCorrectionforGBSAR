function Lo=Loxiangguan(A,B)

if size(A)~=size(B)
    error('A\B 大小不一样');
end

r=size(A,2);
sum_A=sum(sum(A));
sum_B=sum(sum(B));
sum_A_B=sum(sum(A.*B));
sum_A_A=sum(sum(A.*A));
sum_B_B=sum(sum(B.*B));



Lo=(sum_A_B-sum_A*sum_B/r/r)/sqrt((sum_A_A-sum_A*sum_A/r/r)*(sum_B_B-sum_B*sum_B/r/r));
%Lo=sum_A_B/sqrt(sum_A_A*sum_B_B);
