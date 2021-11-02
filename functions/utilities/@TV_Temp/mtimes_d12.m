function res = mtimes(a,b)

if a.adjoint
    res = adjDz(b);
else
    res = b(:,13:end) - b(:,1:end-12);
end

function y = adjDz(x)
mm = size(x, 1);
zzeros = zeros(mm, 12);
y= [zzeros x] - [x zzeros];