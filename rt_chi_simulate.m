function k = rt_chi_simulate(chi, P0, M0, n, asymp)
%rt_CHI_SIMULATE Simulate chi-matrix measurement with multinomial random
%variable generator
%   INPUT:
%   chi - chi-matrix of the process to reconstruct
%   P0 - input density matrices
%       P0(:,:,i) - density matrix of the i-th input state
%   M0 - measurement operators for the state at the output of the gate
%       M0{i}(:,:,j) - j-th operator that corresponds to the i-th
%       measurement, sum(M0{i},3) should be equal to identity matrix
%   n - number of measurements (see description of rt_chi_reconstruct)
%   asymp - if true function will return asymptotic values of n (default: false)
%
%   OUTPUT:
%   k - number of counts (see description of rt_chi_reconstruct)

E = rt_chi2kraus(chi);

mp = size(P0,3);
mm = length(M0);
m = mp*mm;

s = size(P0,1);

if length(n) == 1
    n = ones(1,m)*n;
elseif m ~= length(n)
    error('The length of n should be equal to the number of measurement schemes');
end

if nargin < 5
    asymp = false;
end

k = zeros(s,m);
for i = 1:mp
    for j = 1:mm
        p = zeros(s,1);
        for ii = 1:s
            p(ii) = abs(trace( M0{j}(:,:,ii)*rt_kraus(P0(:,:,i),E) ));
        end
        p = p/sum(p);
        
        ind = (i-1)*mm + j;
        if asymp
            k(:,ind) = n(ind) * p';
        elseif s == 2
            k(1,ind) = binornd(n(ind), p(1));
            k(2,ind) = n(i) - k(1,ind);
        else
            k(:,ind) = mnrnd(n(ind), p)';
        end 
    end
end

end

