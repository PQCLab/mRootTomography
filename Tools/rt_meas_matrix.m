function B = rt_meas_matrix(M)

B = reshape(permute(M,[3,2,1]), size(M,3), []);

end

