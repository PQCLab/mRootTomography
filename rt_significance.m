function [pval, chi2, df, n_observed, n_expected] = rt_significance(dm, clicks, proto, varargin)
%RT_SIGNIFICANCE TODO Vary df

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'dm');
addRequired(p, 'clicks');
addRequired(p, 'proto');
addOptional(p, 'nshots', 'sum');
addParameter(p, 'fromClicks', true);
addParameter(p, 'rank', 'dm');
addParameter(p, 'isProcess', false);
addParameter(p, 'normalizeDM', true);
addParameter(p, 'measDF', 'proto');
parse(p,dm,clicks,proto,varargin{:});
opt = p.Results;

[proto,nshots] = rt_proto_check(proto,opt.nshots,clicks);
[M, n, n_observed] = rt_data_join(proto, nshots, clicks);
n_expected = real(rt_meas_matrix(M)*dm(:)).*n;
if opt.normalizeDM
    n_expected = n_expected / sum(n_expected) * sum(n_observed);
end

if ischar(opt.measDF) && strcmpi(opt.measDF,'proto')
    s = size(dm,1);
    if opt.isProcess
        ss = sqrt(s);
        nPovm = sum(cellfun(@(X) norm(rt_prttrace(sum(X,3),[ss,ss],1)-eye(ss))<1e-5, proto));
    else
        nPovm = sum(cellfun(@(X) norm(sum(X,3)-eye(s))<1e-5, proto));
    end
    df = length(n_observed) - nPovm - (opt.normalizeDM && (nPovm < length(proto)));
else
    df = opt.measDF;
end

if opt.fromClicks
    r = opt.rank;
    if ischar(r) && strcmpi(r,'dm')
        r = rank(dm);
    end
    if opt.isProcess
        nuP = (2*s-r)*r - s;
    else
        nuP = (2*s-r)*r - 1;
    end
    df = df - nuP;
end

chi2 = sum((n_expected-n_observed).^2./n_expected);

if df == 0
    pval = nan;
else
    pval = gammainc(chi2/2,df/2,'upper'); % same as 1-chi2cdf(chi2, df);
    if pval == 0
        pval = vpa(gammainc(chi2/2,df/2,'upper'), 100); % need more precision
    end
end

end
