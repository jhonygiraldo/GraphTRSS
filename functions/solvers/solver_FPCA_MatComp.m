%--------------------------------------------------------------------------
% Fixed Point Continuation Methods based on Approximate SVD (FPCA) for Matrix
% Completion Problem. 
%
% Solves
%           min  rank(X)
%           s.t. X_ij = M_ij , (i,j) \in \Omega
%
% Please refer to One_run.m to get information on how to use this solver.
%
% Reference: "Fixed point and Bregman iterative methods for matrix
%             rank minimization." 
%             Shiqian Ma, Donald Goldfarb, Lifeng Chen, 
%             Technical Report, Columbia University, October 2008
% Available at: 
%             http://www.optimization-online.org/DB_HTML/2008/11/2151.html
%             ftp://ftp.math.ucla.edu/pub/camreport/cam08-78.pdf
% 
% Author: Shiqian Ma 
% Copyright (C) 2008 Shiqian Ma, Columbia University 
% Date: May,15,2008. (Last Modified: July, 06, 2008)
%--------------------------------------------------------------------------
function Out = solver_FPCA_MatComp(m,n,A,b,opts)

% calculate AtMb
Atb = AtMVMx(A,b,m,n); % p = length(b); 
opts.x0 = zeros(m,n); mu = opts.mu; % opts.x0 = opts.tau*Atb;

if max(m,n) > 1000
    [U,sigma] = LinearTimeSVD(Atb,min([m/2,n/2,1000]),min([m/2,n/2,3]),ones(n,1)/n);
    nrm2Atb = sigma(1);
else
    nrm2Atb = norm(Atb);
end


muf = mu; % final value of mu
x = opts.x0; tau = opts.tau;
mu = nrm2Atb*opts.eta; if mu < muf, mu = muf; end

innercount = 0; nrmxxp = inf; count_nrmxxp_increase = 0; extra_rank = 0;
pp = ones(n,1)/n; sn = min([m/2,n/2,round(2*opts.maxr-2)]);
g = get_g(x,m,n,A,Atb);
% main loop
for i = 1:opts.mxitr
    xp = x; nrmxxpp = nrmxxp;
    y = x - tau*g;
    if i == 1
        xp_good = false;
        if max(m,n) < 1000
            [U,S,V] = svd(y); sigma = diag(S);
        else
            [U,sigma] = LinearTimeSVD(y,sn,opts.maxr,pp);
            invsigma = zeros(size(sigma)); indx = find(sigma);
            invsigma(indx) = 1./sigma(indx);
            V = (diag(invsigma)*U'*y)';
        end
    else
        
        sp = s(s>0); mx = max(sp);
        kk = length(find(sp > mx*opts.fastsvd_ratio_leading));
        kk = max(1,min(kk+extra_rank,sn));
        [U,sigma] = LinearTimeSVD(y,sn,kk,pp);
        invsigma = zeros(size(sigma)); indx = find(sigma);
        invsigma(indx) = 1./sigma(indx);
        V = (diag(invsigma)*U'*y)';
    end
    s = sigma;
    ind = find(s > 0);
    Ue = U(:,ind); Ve = V(:,ind); s = s(ind); S = diag(s);

    nu = tau*mu;
    S = max(0,S-nu);
    x = Ue*S*Ve';
    s = diag(S);
    g = get_g(x,m,n,A,Atb);

    nrmxxp = norm(x - xp,'fro');

    if xp_good == true
        if  nrmxxp > nrmxxpp
            count_nrmxxp_increase = count_nrmxxp_increase + 1;
        end
        if count_nrmxxp_increase >= 10
            extra_rank = 1;
            count_nrmxxp_increase = 0;
        else
            extra_rank = 0;
        end
    end
    xp_good = true;
    critx = nrmxxp/max(norm(xp,'fro'),1); innercount = innercount + 1;
    if mu == muf
        opts.maxinneriter = 500;
    end

    if (critx < opts.xtol) || (innercount >= opts.maxinneriter)
        innercount = 0; xp_good = false;
        % stop if reached muf
        if mu == muf
            Out.x = x; Out.iter = i;
            return
        end

        mu = opts.eta*mu;
        mu = max(mu,muf);
    end

end

% did not converge within opts.mxitr
Out.x = x; Out.iter = i;

%--------------------------------------------------------------------------
% SUBFUNCTION FOR CALCULATING g
%--------------------------------------------------------------------------
    function y = AMVMx(A,x)
        % y = A*x
        y = x(A);
    end

    function y = AtMVMx(A,b,m,n)
        % y = A'*b
        y = zeros(m*n,1);
        y(A) = b;
        y = reshape(y,m,n);
    end


    function g = get_g(x,m,n,A,Atb)
        Ax = AMVMx(A,reshape(x,m*n,1));
        g = AtMVMx(A,Ax,m,n)-Atb;
    end
end