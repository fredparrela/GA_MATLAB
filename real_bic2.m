function [dic] = real_bic2(numParam,count,obs,dic,stado_pai,stado_variavel,names)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% tamanho=0;
% N_obs=CCC;
 bic1=zeros(size(numParam,1),1);
% serie=[];
% var=[];

for n=1:count
    % bic1(n)=0;
    %     names = string.empty(0,0);
    para=(numParam{n});
    %    para=para+exp(-700);


    %     if(size(para,2)==3)
    % %         stado_pai=1;
    % %         stado_variavel(n)=size(para,1);
    %
    % %      names=(para.Properties.VariableNames{1});
    %     else
    %         aux=string(para.Properties.VariableNames);
    %         aux=aux(1:end-2);
    %         aux3 = '';
    %         for kk=1:size(aux,2)
    %             aux3=aux(kk)+aux3;
    %         end
    %         names=aux3;
    % %         stado_pai=height(unique(para(:,1:end-3),'rows'));
    % %         stado_variavel(n)=height(unique(para(:,end-2)));
    %     end

    %     tamanho=tamanho+(stado_variavel(n)-1)*stado_pai;




    for q=1:stado_pai(n)


        if(size(para,2)==3)
            aux=para.GroupCount;
            paraLn=sum(aux,"all");
        else

            if(q==1)
                aux=para.GroupCount(1:stado_variavel(n));
            else
                aux=para.GroupCount((q-1)*stado_variavel(n) +1 : (q-1)*stado_variavel(n) + stado_variavel(n));

            end
            %           aux=aux+exp(-700);
            paraLn=sum(aux);
        end
        if((paraLn)==0)
            aux=aux+1;
            paraLn=sum(aux);
        end
    
        for r=1:stado_variavel(n)

            % if((paraLn)==0)

                 % bic1(n)=bic1(n)-0.1*log(obs);
                %                 serie=[serie 0];
                %                 var=[var n];
                 % 1+1;
            % else
                if(aux(r)==0)
                    %                     serie=[serie 0];
                    %                     var=[var n];
                    bic1(n)=bic1(n)+0;
                else
                    %                 bic1(n)=bic1(n)+ log(( aux(r)/paraLn ));
                    bic1(n)=bic1(n)+ (aux(r))*log(( aux(r)/paraLn ));
                    %                     serie=[serie (aux(r))*log(( aux(r)/paraLn ))];
                    %                     var=[var n];

                end

            % end

        end




    end

   
    %       serie=[serie bic1(n)];

end
 dic(names)=bic1';
% bic1(n)=bic1(n)-0.5*log(obs)*tamanho ;

end

