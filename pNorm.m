function [norm_result] = pNorm(A_matrix, norm_value)
% Calculates the (p) norm, if inf use a very large number

sum = 0;
for ind_i = 1:size(A_matrix, 1)
    for ind_k = 1:size(A_matrix, 2)
        sum = sum + abs(A_matrix(ind_i, ind_k))^norm_value;
    end
end
norm_result = sum^(1 / norm_value);
end
