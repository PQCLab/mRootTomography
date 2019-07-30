function E = rt_kraus_matrices(type, param, q_num, q_size)
%KRAUS_MATRICES Returns Kraus matries with specified parameters

E = [];
if strcmpi(type, 'dephasing')
    if length(param) == 2
        param = param(1) / param(2);
    end
    gamma = 1 - exp(-2*param);
    E(:,:,1) = [1 0; 0 sqrt(1-gamma)];
    E(:,:,2) = [0 0; 0 sqrt(gamma)];
elseif strcmpi(type, 'amplitude')
    if length(param) == 2
        param = param(1) / param(2);
    end
    gamma = 1 - exp(-param);
    E(:,:,1) = [1 0; 0 sqrt(1-gamma)];
    E(:,:,2) = [0 sqrt(gamma); 0 0];
elseif strcmpi(type, 'relaxation')
    if length(param) < 3
        param = [1 param];
    end
    T = param(1);
    T1 = param(2);
    T2 = param(3);
    T2_pure = 2*T1*T2 / (2*T1 - T2);
    
    Ea = kraus_matrices('amplitude', [T,T1]);
    Ep = kraus_matrices('dephasing', [T,T2_pure]);
    
    e = [reshape(Ea(:,:,1)*Ep(:,:,1),4,1),...
        reshape(Ea(:,:,1)*Ep(:,:,2),4,1),...
        reshape(Ea(:,:,2)*Ep(:,:,1),4,1),...
        reshape(Ea(:,:,2)*Ep(:,:,2),4,1)]/sqrt(2);
    
    E = chi2kraus(e*e');
    
elseif sum(strcmpi(type, {'X', 'Y', 'Z'}))
    E(:,:,1) = sqrt(1-param) * eye(2);
    E(:,:,2) = sqrt(param) * pauli(type);
elseif strcmpi(type, 'depolarizing')
    E(:,:,1) = sqrt(1-3*param/4) * eye(2);
    E(:,:,2) = sqrt(param/4) * pauli('x');
    E(:,:,3) = sqrt(param/4) * pauli('y');
    E(:,:,4) = sqrt(param/4) * pauli('z');
elseif strcmpi(type, 'random')
    if length(param) == 1
        s = 2;
        r = param;
    else
        s = param(1);
        r = param(2);
    end
    U = rt_randunitary(s*r);
    E = zeros(s,s,r);
    for i = 1:r
        E(:,:,i) = U((i-1)*s+(1:s), 1:s);
    end
elseif strcmpi(type, 'chi')
    E = rt_chi2kraus(param);
else
    error('Unknown type');
end

if nargin == 3
    q_size = q_num;
    E = chi2kraus(kraus2chi(kronpower(E, q_size)));
elseif nargin == 4
    Eq = E;
    r = size(Eq,3);
    E = zeros(2^q_size, 2^q_size, r);
    for k = 1:r
        E(:,:,k) = kronm(eye(2^(q_num-1)), Eq(:,:,k), eye(2^(q_size - q_num)));
    end
end

end
