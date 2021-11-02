function x_recon = solver_NNI(SampleMatrix, SampledData, Position)
%
% natural neighbor interpolation
% 
%  Inputs:
%
%    SampleMatrix -- the sampling operator
%   
%    SampledData  -- the sampled data
%
%    Position     -- the positions of all the vertices of the graph
%
%  Output:
%
%      x_recon    -- the reconstructed graph signal
% 
%  Written by Kai Qiu (q1987k@163.com)
%  date 6/20/2017


[N,T] = size(SampledData);
x_recon = SampledData;
xi=[];yi=[];ti=[];zi=[];
for TimeInd=1:T
    temp1=SampleMatrix(:,TimeInd);
    temp2=SampledData(:,TimeInd);
    x1=find(temp1~=0);
    z1=temp2(x1);
    xi=[xi;Position(x1,1)];
    yi=[yi;Position(x1,2)];
    ti=[ti;TimeInd*ones(length(x1),1)];
    zi=[zi;z1];
end

F=scatteredInterpolant(xi,yi,ti,zi,'natural'); % nearest linear natural
for TimeInd=1:T
    temp1=SampleMatrix(:,TimeInd);
    x2=find(temp1==0);
    x_recon(x2,TimeInd)=F(Position(x2,1),Position(x2,2),TimeInd*ones(length(x2),1));
end
