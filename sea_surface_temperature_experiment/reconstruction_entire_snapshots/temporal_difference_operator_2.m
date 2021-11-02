function D_h = temporal_difference_operator_2(signals_t)
% This function returns the temporal difference operator with 2 steps.
D_h = zeros(signals_t,signals_t-2);
D_h(1,1) = -1;
D_h(2,1) = 1;
D_h(2,2) = -1;
D_h(end,end) = -1;
D_h(end-1,end) = 1;
D_h(end-1,end-1) = -1;
for i=3:signals_t-2
    D_h(i,i-1) = 1;
    D_h(i,i) = -1;
    D_h(i,i-2) = -1;
end