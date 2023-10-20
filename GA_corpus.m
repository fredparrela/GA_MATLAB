close all
clearvars
%% Import data

% profile on

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

pop_ini_size=400;

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
best_fit={};
pop_ini={};
%%
parfor i=1:10
    [best_fit{i},pop_ini{i}] =GA(df,n_pais,PC,PM,PO,pop_ini_size,n,dic,summ,prime_list,obs,index,size_bit);
end
[fit_truth,dic] = score_ot_full(dag_truth,df,obs,n,dic,summ,prime_list);
%%
for i=1:10
    % figure
    % imagesc(pop_ini{i})
    % figure
    % plot(best_fit{i})
    dag=bit_to_dag( pop_ini{i}(1,:),n,index);
    digraph_best=digraph(dag,topological_order);
    figure
    plot(digraph_best)
    score_diff(i)=best_fit{i}(end)-fit_truth;
    title(score_diff(i))
    [BSF(i)] = Com_resul(dag,dag_truth,n);
    %display(BSF)
end
figure
bar(BSF)
figure
bar(score_diff)
% profile off
% profile viewer
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



