function [dag] = bit_to_dag(bit_represent,n,index)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
F1=zeros(n);


% rows=[];
% coluns=[];
% for i=1:n-1
% 
%     count=ini:ini+(n-i-1);
%     F1(i,i+1:n)=bit_represent(count);
%     ini=count(end)+1;
%     % rows=[rows i];
%     % coluns=[coluns i+1:n];
% 
% 
% 
% end

for i=1:size(bit_represent,2)

F1(index(i,1),index(i,2))=bit_represent(i);

end


%%

dag=1*((F1>0) +(F1'<0));
end