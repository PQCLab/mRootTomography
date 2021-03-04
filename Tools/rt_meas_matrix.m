function B = rt_meas_matrix(P)

B = reshape(permute(P,[3,2,1]), size(P,3), []);

end

