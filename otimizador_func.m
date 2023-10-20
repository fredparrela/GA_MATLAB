function [dag,Fit_1,dic,bit_represent] = otimizador_func(bit_represent,DATA_2,obs,n,dic,summ,primos,index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dag=bit_to_dag(bit_represent,n,index);

[Fit_1,dic]=score_ot_full(dag,DATA_2,obs,n,dic,summ,primos);


[bit_represent2] = mutation_bit_2(bit_represent,n,index);
[dag2] = bit_to_dag(bit_represent2,n,index);

[Fit_2,dic,~]=score_ot_full(dag2,DATA_2,obs,n,dic,summ,primos);
busca=0;
k=0;
while( and(k<n,1))

    if(Fit_2>Fit_1 )

        dag=dag2;
        busca=1+busca;
        Fit_1=Fit_2;
        bit_represent=bit_represent2;

    else

        [bit_represent2] = mutation_bit_2(bit_represent,n,index);
        [dag2] = bit_to_dag(bit_represent2,n,index);
        [Fit_2,dic,~]=score_ot_full(dag2,DATA_2,obs,n,dic,summ,primos);

    end

    k=k+1;
end


end