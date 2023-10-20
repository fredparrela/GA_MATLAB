function [best_fit,pop_ini] =GA(df,n_pais,PC,PM,PO,pop_ini_size,n,dic,summ,prime_list,obs,index,size_bit)
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
end