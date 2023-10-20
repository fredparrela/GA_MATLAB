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
