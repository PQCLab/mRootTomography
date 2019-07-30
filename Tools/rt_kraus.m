function rho = rt_kraus(rho, E)
%KRAUS performs kraus transform

rho_in = rho;
if size(rho_in, 2) == 1
    rho_in = rho_in*rho_in';
end

rho = 0;
for k = 1:size(E,3)
    rho = rho + E(:,:,k)*rho_in*E(:,:,k)';
end

end

