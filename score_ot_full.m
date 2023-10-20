function [bic11,dic,tamanho] = score_ot_full(F1,DATA,obs,n,dic,summ,primos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pred={};
  names = zeros(1,size(F1,1));
 %names=string([]);

fn = fieldnames(summ);
stado_variavel=zeros(n,1);stado_pai=zeros(n,1);
for node=1:n
    pais=find(F1(:,node)==1);
    pred{node}=pais;
    %names(node)=append(DATA.Properties.VariableNames{node},DATA.Properties.VariableNames{pais});
    % [names(node)] = name_order(DATA,pais,node);
     names(node)=name_order_2(pais,node,n,primos);
    stado_variavel(node)=size(summ.(fn{node}).Categories,1);
    if(isempty(pais))
        stado_pai(node)=1;
    else
        aux=1;
        for k=1:size(pais,1)
            aux=aux*size(summ.(fn{pais(k)}).Categories,1);
        end
        stado_pai(node)=aux;
    end
end

tamanho=sum((stado_variavel-1).*stado_pai);

tf = isKey(dic,names);
colu=find((~tf));
% pred_1=pred((~tf));
if(isempty(colu))

else
    [test,count]=Loglikehood_g8(pred,DATA,colu,fn);
    [dic] = real_bic2(test,count-1,obs,dic,stado_pai(~tf),stado_variavel(~tf),names(~tf));
end

% if(tamanho>509)
%     bic11=sum(dic(names))-(tamanho-509)*log(obs)*tamanho ;
% else
% 
%     bic11=sum(dic(names))-0.5*log(obs)*tamanho ;
% end

 bic11=sum(dic(names))-0.5*log(obs)*tamanho ;
%[bic11] = real_Bdeu(test,count-1,obs);
% "SHUNT HREKG HRBP ERRLOWOUTPUT MINVOL INTUBATION" Normal
% "SHUNT MINVOL HREKG HRBP ERRLOWOUTPUT INTUBATION"
% load('test_p.mat')
% load('test_n.mat')

end