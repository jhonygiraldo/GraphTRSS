function [sampling_operator, sampled_indexes] = sampling_chen_algorithm(V_K,m,N)
%This function computes the sampling operator with Chen's algorithm from
%the paper: Discrete Signal Processing on Graphs, Sampling Theory.
%By: Jhony H. Giraldo
possible_indexes = [1:N];
sampled_indexes = [];
sampling_operator = [];
while(length(sampled_indexes) < m)
    length(sampled_indexes)
    smallest_singular_values = zeros(length(possible_indexes),1);
    for(i=1:length(smallest_singular_values))
        M = zeros(1,N);
        M(possible_indexes(i)) = 1;
        sampling_operator_temp = [sampling_operator;M];
        obj_func = sampling_operator_temp*V_K;
        [U,S,V] = svd(obj_func);
        if(size(S,1) == 1)
            smallest_singular_values(i) = S(1,1);
        else
            smallest_singular_values(i) = min(diag(S));
        end
    end
    index_max = find(smallest_singular_values == max(smallest_singular_values));
    sampled_indexes = [sampled_indexes;possible_indexes(index_max(1))];
    possible_indexes = setdiff(possible_indexes,possible_indexes(index_max(1)));
    %%
    M = zeros(1,N);
    M(sampled_indexes(end)) = 1;
    sampling_operator = [sampling_operator;M];
end