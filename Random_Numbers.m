function RN = Random_Numbers(lower_bound, upper_bound, n_row, n_column)
temp = rand(n_row, n_column);
RN = (upper_bound - lower_bound) .* temp + lower_bound;
end