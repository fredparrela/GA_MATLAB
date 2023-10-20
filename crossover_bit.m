function [F1,F2] = crossover_bit(A,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%  N1=size(A,1);

%  colum=randperm(N1,2);
%
F2=zeros(size(A));
F1=zeros(size(A));
for k=1:size(A,1)
    %
    AA=A(k,:);
    BB=B(k,:);
    %%

    point=randperm(numel(AA),2);
    point=sort(point);
    %
%     point=round((n*n)/2);
%     aux1=reshape([AA(1:point(1)) BB(point(1)+1:end)],size(AA));
%     aux2=reshape([BB(1:point(1)) AA(point(1)+1:end)],size(AA));


    aux1=reshape([AA(1:point(1)) BB(point(1)+1:point(2)) AA(point(2)+1:end)],size(AA));


    aux2=reshape([BB(1:point(1)) AA(point(1)+1:point(2)) BB(point(2)+1:end)],size(AA));

    %     Cr1=AA(mask==1);
    %     Cr2=BB(mask==1);
    %
    %     aux=or(Cr1,Cr2);
    %
    %
    %
    %     AA(mask==1)=(aux.*rand(sum(mask==1,'all'),1)  >0.5);
    %     BB(mask==1)=(aux.*rand(sum(mask==1,'all'),1)  >0.5);
    %
    %%


    F2(k,:)=aux1;
    F1(k,:)=aux2;

end

% F2= (A.*B);
% F1= or(A,B);

% A B
% 1 1 = 0
% 1 0 = 1
% 0 1 = 1
% 0 0 = 0;
end