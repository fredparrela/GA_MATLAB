function [bit_represent] = dag_to_bit(F1,n)
%Summary of this function goes here
%   Detailed explanation goes here
bit_represent=zeros(1,n*(n-1)*.5);
ini=1;
F=triu(F1,1)+ (-1*tril(F1,-1))';

for i=1:n-1

    count=ini:ini+(n-i-1);
    bit_represent(count)=F(i,i+1:n);
    ini=count(end)+1;

end
end
