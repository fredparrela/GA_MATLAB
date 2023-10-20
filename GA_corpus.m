close all
clearvars
%% Import data

profile on

filename='ALARM';
filename_data=strcat(filename,'_DATA.csv');
obs=50000;
df=import_data(filename_data);
sort_rows=randperm(size(df,1),obs);
n=(size(df,2));


top_sort=randperm(n,n);
df=df(sort_rows,top_sort);

topological_order=df.Properties.VariableNames;


filename_dag=strcat('DAGtrue_',filename,'.csv');
dag_truth=read_dag_true(filename_dag,topological_order);

di_graph=digraph(dag_truth,topological_order);

figure
plot(di_graph)
node_names=topological_order;
dic=dictionary([],[]);
prime_list=(primes(2000));
%% GA parameters
summ=summary(df);

pop_ini_size=350;

size_bit=((n-1)*n*0.5);
PC=0.5;
PM=1;
PO=0.9;
n_pais=10;
%%
mask=triu(ones(n),1);
[rows,cols]=find(mask==1);
index =[rows cols];
index=sortrows(index);
pais=zeros(n_pais,1);
%%

pd = makedist('Multinomial','Probabilities',[1/size_bit (size_bit-2)/size_bit 1/size_bit]);
[pop_ini] = Pop_ini_bit(pop_ini_size,pd,n);
[dag_pop, fit]=Pop_ini_dag (pop_ini,df,pop_ini_size,...
    n,index,obs,dic,summ,prime_list);

[fit, index_fit]=sort(fit,'descend');
dag_pop=dag_pop(:,:,index_fit);
pop_ini=pop_ini(index_fit,:);

%% selecao de pais
geracao=1;
stable_count = 0;
stability_threshold = 500;
max_iterations = 6000;
prev_best_cost=0;
best_fit(geracao)=fit(1);
%%
while(1)
    filhos=[];
    selecao_torneio=randperm(pop_ini_size,2*n_pais);
    [pais]=sort(selecao_torneio);
    pais=pais((1:n_pais));
    fit_filhos=zeros(n_pais,1);

    for i=1:(n_pais/2)
        if(PC>rand)
            [f1,f2]=crossover_bit( pop_ini(pais(i),:) ,pop_ini(pais(i+1),:)  );

            [dag_f1] = bit_to_dag(f1,n,index);
            Gi= digraph(dag_f1);
            if(~isdag(Gi))
            [dag_f1] = repair_dag(dag_f1);
            end
            

            [dag_f2] = bit_to_dag(f2,n,index);
            Gi= digraph(dag_f2);
            if(~isdag(Gi))
            [dag_f2] = repair_dag(dag_f2);
            end

            [f1] = dag_to_bit(dag_f1,n);
            [f2] = dag_to_bit(dag_f2,n);

            filhos=[filhos; f1;f2];
        else
            filhos=[filhos; pop_ini(pais(i),:); pop_ini(pais(i),:)];
        end

    end

    for i=1:n_pais
        if(PM>rand)
            filhos(i,:)=mutation_bit_2(filhos(i,:),n,index);
        end
    end
    if(PO>rand)
        idx=randperm(n_pais,1);
        [~,Fit_1,dic,bit_represent] = otimizador_func(filhos(idx,:),df,...
            obs,n,dic,summ,prime_list,index);
        fit_filhos(idx)=Fit_1;
        filhos(idx,:)=bit_represent;
    end

    for i=1:n_pais
        dag=bit_to_dag(filhos(i,:),n,index);
        %if(fit_filhos(i)==1)
        [fit_filhos(i),dic] = score_ot_full(dag,df,obs,n,dic,summ,prime_list);
        %end
    end

    fit_final_bin=[fit ;fit_filhos];
    [~,I]=sort(fit_final_bin,'descend');
    pop_aux=[pop_ini ;filhos];

    pop_ini=pop_aux(I(1:pop_ini_size),:);
    fit=fit_final_bin(I(1:pop_ini_size));


    [fit,I]=sort(fit,'descend');
    pop_ini=pop_ini(I,:);

    geracao=geracao+1;

    best_fit(geracao)=fit(1);
    prev_best_cost=best_fit(geracao-1);
    
    
    

    % Check if the current best_fit is equal to the previous best_cost
    if best_fit(geracao) == prev_best_cost
        stable_count = stable_count + 1;
        % disp(['stable_count: ' num2str(stable_count)]);
    else
        stable_count = 0;
    end

    % Check for stability
    if stable_count >= stability_threshold
        disp(['Stopping criterion met: Cost has been stable for ' num2str(stability_threshold) ' generations.']);
        break;
    end

    % Check if the maximum number of iterations has been reached
    if geracao >= max_iterations
        disp('Maximum number of iterations reached.');
        break;
    end


    
end
%%
[fit_truth,dic] = score_ot_full(dag_truth,df,obs,n,dic,summ,prime_list);

figure
imagesc(pop_ini)
figure
plot(best_fit)
dag=bit_to_dag( pop_ini(1,:),n,index);
digraph_best=digraph(dag,topological_order);
figure
plot(digraph_best)
title(best_fit(end)-fit_truth)

[BSF] = Com_resul(dag,dag_truth,n);
display(BSF)

profile off
profile viewer
%%
% function my_dict = my_dic_init(topological_order)
%     % Generate a dictionary (map) that maps column names to their indices.
%     % Input:
%     % - data: Data table or dataset
%     % Output:
%     % - my_dict: Dictionary mapping column names to indices
%
%     my_dict = dictionary(string([]),[]);
%
%     for index = 1:length(topological_order)
%         my_dict(topological_order{index}) = index;
%     end
% end
function [dag_pop, fit]=Pop_ini_dag (pop_ini,df,pop_ini_size,...
    n,index,obs,dic,summ,prime_list)
dag_pop=zeros(n,n,pop_ini_size);
fit=zeros(n,1);
for i=1:pop_ini_size
    bit_represent=pop_ini(i,:);
    [dag] = bit_to_dag(bit_represent,n,index);
    dag_pop(:,:,i)=dag;
    Gi= digraph(dag);
    if(~isdag(Gi))
        [dag] = repair_dag(dag);
        [bit_represent] = dag_to_bit(dag,n);
        pop_ini(i,:)=bit_represent;
    end
    [bic11,dic] = score_ot_full(dag,df,obs,n,dic,summ,prime_list);
    fit(i)=bic11;

end

end

function [filhos_bin] = repair_dag(filhos_bin)
%UNTITLED2 Summary of this function goes here

for pp=1:size(filhos_bin,3)
    aux=filhos_bin(:,:,pp);
    Gi= digraph(aux);

    while(~isdag(Gi))
        [~,edgecycles] = allcycles(Gi);
        k=randperm(size(edgecycles,1),1);
        Gi=rmedge(Gi,edgecycles{k});
        if(isdag(Gi))
            filhos_bin(:,:,pp)=full( adjacency(Gi));
            break

        end


    end

end
end

function table=import_data(filename)
% Create a new import options object
import_options = detectImportOptions(filename);
% Set the VariableTypes property to 'categorical'
n=size(import_options.VariableTypes,2);
for i=1:n
    import_options.VariableTypes{i} = 'categorical';
end
% Read the CSV file into a table
table = readtable(filename, import_options);

%%
end

function adjacency_matrix=read_dag_true(filename,topological_order)
%% Read the CSV file into a table
df = readtable(filename);

s=df(:,2).Variables;
t=df(:,4).Variables;
G= digraph(s,t,[],topological_order);
A = adjacency(G,'weighted');
adjacency_matrix= full(A);


end

% function group_counts = state_count(data, adjacency_matrix, node, node_names)
%
% % Compute the states conditionally according to the DAG represented by the adjacency matrix.
% %
% % Args:
% %
% % adjacency_matrix: A adjacenty matrix.
% % names:The coluns, nodes name.
% % data: Data frame ( All types must be converted into categorical)
% % node: The colum/ node which the states will be computed
% %
% % Returns:
% %  data frame containing the states conditionally counted according to the adjacency matrix
% %
%
%     parents = adjacency_matrix(:, node) == 1;
%     selected_node_names = node_names(parents);
%     group_filtered = [selected_node_names, node_names(node)];
%
%     %group_counts = grpstats(data, group_filtered, 'size', 'DataVars', false);
%     group_counts=groupcounts(data,group_filtered,"IncludeEmptyGroups",true);
%     group_counts=group_counts(:,1:end-1);
% end

% function [bic_score,my_dict_global]
% my_dict = my_dic_init(topological_order);
% dic=dictionary([],[]);
% profile on
% [bic_score,dic]=BIC_score(df, dag_truth, node_names,prime_list,dic,my_dict);
% profile off
% profile viewer
%
% function [bic_score,my_dict_global] = BIC_score(data, adjacency_matrix, node_names,prime_list,my_dict_global,my_dict)
%     % Compute the BIC score of DAG given the data, adjacency matrix, and node names.
%     % Input:
%     % - data: Data matrix (all types must be converted into categorical)
%     % - adjacency_matrix: An adjacency matrix.
%     % - node_names: The columns, nodes' names.
%     % Output:
%     % - BIC score of the DAG
%
%     n = size(adjacency_matrix, 1);
%     obs = size(data, 1);
%     Bn_size = 0;
%     acc = 0;
%
%     for node = 1:n
%         bic = 0;
%         states = state_count(data, adjacency_matrix, node, node_names);
%         teste = states.Properties.VariableNames(1:end-1);
%         key = 1;
%
%         if length(teste) == 1
%             key = prime_list(my_dict(teste{1}) + 2 * n);
%         else
%             aux = teste(1:end-1);
%             for i = 1:length(aux)
%                 key = key * prime_list(my_dict(aux{i}));
%             end
%             key = key * prime_list(my_dict(teste{end}) + 2 * n);
%         end
%
%         A = length(unique(states{:, end - 1}));
%         B = size(states, 1);
%         Bn_size = (B / A) * (A - 1) + Bn_size;
%
%         if ~isKey(my_dict_global, key)
%             count = 0;
%             count2 = 0;
%
%             for i = 1:(B / A)
%                 for z = 1:A
%                     if sum(states{count2 + 1:count2 + A, end}) == 0
%                         states{count2 + 1:count2 + A, end} = 1;
%                         bic = bic + states{count + 1, end} * log(states{count + 1, end} / sum(states{count2 + 1:count2 + A, end}));
%                     else
%                         if states{count + 1, end} == 0
%                             bic = bic + 0;
%                         else
%                             bic = bic + states{count + 1, end} * log(states{count + 1, end} / sum(states{count2 + 1:count2 + A, end}));
%                         end
%                     end
%                     count = count + 1;
%                 end
%                 count2 = count2 + A;
%             end
%
%             my_dict_global(key) = bic;
%         end
%
%         acc = acc + my_dict_global(key);
%     end
%
%     bic_score = acc - 0.5 * log(obs) * Bn_size;
% end

function [pop_ini] = Pop_ini_bit(pop_ini_size,pd,n)


G_lib=(n*(n-1)*0.5);
pop_ini=zeros(pop_ini_size,G_lib);
bit_represent=zeros(G_lib,1);
for k=1:pop_ini_size
    ini=1;
    for i =1:n-1
        count=ini:ini+(n-i-1);
        bit_represent(count)= random(pd,1, n-i);
        ini=count(end)+1;
    end
    bit_represent=bit_represent-2;
    pop_ini(k,:)=bit_represent;
end

end

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