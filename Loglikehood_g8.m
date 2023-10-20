function [state,count] = Loglikehood_g8(G,DATA,Colu,fn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
count=1;
state={};
% tiny = exp(-700);


for qu =1:size(Colu,2)
    i=Colu(qu);
    if(isempty(G{i}))

        [GC]= groupcounts(DATA(:,i),fn{i},"IncludeEmptyGroups",true);
        state{count,1}= GC;
        % state{count,2}=1;
        count=count+1;



    else
        pre=G{i};
        pre=pre';
        %         V={V1,V2,V3};
        %         cat_parent={};
        Names={};

        for p=1:size(pre,2)


            Names{p}=(fn{pre(p)});
            %              cat_parent{p}=summ.(fn{pre(p)}).Categories;

        end


        %         cat_parent{p+1}=summ.(fn{i}).Categories;
        Names{p+1}=(fn{i});
        %         [index,~] = cartesianProd(cat_parent);

        %         estado=[];
        % %         tic
        %         final=zeros(size(index,1),1);
        %          tic
        %         for k=1:size(index,1)
        %
        %             for kk=1:size(pre,2)+1
        %                 estado=[estado  string(prd{k,kk}{1})];
        %             end
        %
        %
        %              final(k)=sum(prod(DATA(:,[i pre]).Variables==categorical(estado),2));
        % %               final(k)=sum(ismember(DATA(:,[pre i]).Variables,categorical(estado),'rows'));
        % % %             v=DATA(:,[pre i]).Variables;
        % % %             v2=categorical(estado);
        % %
        % %             out = sum(ismember(DATA(:,[pre i]).Variables,categorical(estado),'rows'));
        %
        %
        %             estado=[];
        %
        %         end
        %                  toc

        [GC]= groupcounts(DATA(:,[i pre]),Names,"IncludeEmptyGroups",true);


        final=GC;
        state{count,1}= (final);

        count=count+1;


    end
end
end
