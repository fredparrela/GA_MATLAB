function [num] = name_order_2(pais,node,n,primos)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%pais=sort(pais);


try     
    num=prod([ primos(n+node) primos(pais')]);
catch
    warning('Problem using function.  Assigning a value of 0.');
    num = n;

end
% n = max(seq); % determine the number of possible values for each digit
% base_n = n.^(0:length(seq)-1); % compute the base-n values for each digit
% num = sum(seq.*base_n); % compute the final integer by summing the products of each digit and its corresponding base-n value

%num= ((num2str(seq)));

% teste2=[node pais'];
end