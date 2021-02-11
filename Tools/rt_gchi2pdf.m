function [p, x] = rt_gchi2pdf(x, d)
%GCHI2PDF TODO

x0 = x;
tol = 1e-4;
nbins = 1e4;

xmax = chi2inv(1-tol, 1) * sum(d);
dx = xmax / nbins;
x = (0:nbins-21) * dx;
nu = length(d);

if nu < 3 % Exact solution for 1 or 2 degrees of freedom
    if ~isempty(x0)
        x = x0;
    end
    if nu == 1
        p = chi2pdf(x / d, 1) / d;
    elseif nu == 2
        a = d(1);
        b = d(2);
        p = 1 / (2 * sqrt(a * b)) * exp(-(a + b) * x / (4 * a * b)) .* besseli(0, (a - b) * x / (4 * a * b));
    end
    return;
end

% Numerical solution for arbitrary degrees of freedom
du = 2 * pi /xmax;
u = (0:nbins-1) * du;
phi = prod(1 ./ sqrt(1 - 2*1j*bsxfun(@times, u, d(:))), 1);
p = (du / pi) * real(fft(phi));
p = p(1:nbins-20) - min(p);

if ~isempty(x0)
    p = interp1(x, p, x0, 'linear', 'extrap');
end

end

