function [Sopt,Sindex]=sampling_tsitsvero_algorithm(U,M,BW)

% Matlab code implementing the sampling technique of the paper: "Signals on Graphs:
% Uncertainty Principle and Sampling" % By: Alejandro Parada-Mayorga

% BW: bandwidth of the signal
% M: Number of samples
% U: basis matrix

U_tilde=( U(:,1:BW) )';
Sindex=[];

while max(size(Sindex))<M
    
    K=min(max(size(Sindex)),BW);    
    V=setdiff([1:size(U,1)],Sindex);
    
    for r=1:1:max(size(V))
        s = svd(U_tilde(:,union(Sindex,V(r))));
        fcost(r)=sum(1./s(1:K).^2); clear s
    end
    [~,auxindex]=min(fcost); clear fcost
    Sindex = [Sindex,V(auxindex(1))]
    %Sindex=union(Sindex,V(auxindex(1)))
     
end    

Sopt=zeros(size(U,1),1);
Sopt(Sindex)=1;


end