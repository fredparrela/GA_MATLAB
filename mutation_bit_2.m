function [A] = mutation_bit_2(A,n,index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
flag=1;
for k=1:size(A,1)
    %% Mutation 1
    while(flag)
        aux=A(k,:);
        coin=rand;
        if(coin<1/3)
            set=find(aux==0);
        elseif(coin<2/3)
            set=find(aux==1);
        elseif(coin>2/3)
            set=find(aux==-1);
        end
        %AA=size(set,2);
        if(rand>0.9)
            if(~isempty(set))
                if(coin<1/3)
                    %mu=floor(n/3);
                    mu=3;
                else
                    mu=randperm(length(set),1);
                end
                bit3=set(randperm(length(set),mu));

            else
                bit3=randperm(size(aux,2),1);
            end
        else

            if(isempty(set))
                bit3=randperm(size(aux,2),1);
            else
                bit3=set(randperm(size(set,2),1));
            end
        end
        %%
        if(size(bit3,2)>50)
        bit3=bit3( randperm(size(bit3,2),50)   );
        end

        for i=1:size(bit3,2)
            bit=bit3(i);
            if(aux(bit)==0)
                old=0;
                if(rand>0.5)
                    aux(bit)=1;

                else
                    aux(bit)=-1;
                end

            elseif(aux(bit)==1)
                old=1;

                if(rand>0.5)
                    aux(bit)=0;
                else
                    aux(bit)=-1;
                end

            else
                old=-1;
                if(rand>0.5)
                    aux(bit)=0;
                else
                    aux(bit)=1;
                end

            end

            [dag] = bit_to_dag(aux,n,index);
            Gi= digraph(dag);
            flag=~isdag(Gi);
            if(flag)
                aux(bit)=old;
            end

        end

        %%

    end
    A(k,:)=aux;
    flag=1;
    count=1;
    bit_1=[];
    %% Mutation 2


end


end