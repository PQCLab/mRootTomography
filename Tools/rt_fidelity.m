function f = rt_fidelity(a,b)
%FIDELITY Calculates the fidelity of quantum states
%   Calculates the fidelity of quantum states

if size(a,2) > 1 && size(b,2) > 1
    % f = real(trace(sqrtm(sqrtm(a)*b*sqrtm(a))))^2;
    a = a / trace(a);
    b = b / trace(b);
    [u,d] = eig(a);
    sd = sqrt(d);
    A = u*sd*u' * b * u*sd*u';
    f = real(sum(sqrt(eig(A)))^2);
elseif size(a,2) > 1
    a = a / trace(a);
    b = b / sqrt(b'*b);
    f = abs(b'*a*b);
elseif size(b,2) > 1
    a = a / sqrt(a'*a);
    b = b / trace(b);
    f = abs(a'*b*a);
else
    f = abs(a'*b)^2;
end

end

