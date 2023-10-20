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