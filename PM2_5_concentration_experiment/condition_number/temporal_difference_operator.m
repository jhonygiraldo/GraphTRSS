function D_h = temporal_difference_operator(signals_t)
% This function returns the temporal difference operator.
D_h = sparse(signals_t,signals_t-1);
D_h(1,1) = -1;
D_h(end,end) = 1;
for i=2:signals_t-1
    D_h(i,i-1) = 1;
    D_h(i,i) = -1;
end