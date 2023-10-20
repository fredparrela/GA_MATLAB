function [BSF] = Com_resul(D,F1,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
G2= digraph(F1);
a=numedges(G2);
V=numnodes(G2);
ii=( (V*(V-1) )*.5) -a;
target=reshape(F1,[n*n,1]);




%%
for k=1:size(D,3)
    y=reshape(D(:,:,k),[n*n,1]);
    C = confusionmat(target',y');

%     BSF(k)= 0.5*((C(2,2)/a) + ( ( ii-C(1,2) ) /ii) - (C(1,2) /ii) - (C(2,1)/a)); 
    P=C(2,2)/(C(2,2)+C(2,1) );
    Re=C(2,2)/(C(2,2)+C(1,2) );
    BSF(k)=2*(P*Re/(Re+P) );
    %     figure
%     plotconfusion(target',y')

end
%%
% [~,index]=sort(BSF,'descend');


%teste=D(I,I,1)

% figure
% plot(G2)
% title('Original')
% 
% G1= digraph(D(:,:,index(1)),ASIADATA.Properties.VariableNames);
% 
% figure
% plot(G1)
% title('Aprendido')
end