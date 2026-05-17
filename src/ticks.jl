function round_to_next_N(x; N=100)
    return ceil(x / N) * N
end

function get_integer_ticks(x_lim; N_spacing=100)
    round_to_next_N(x_l[1], N=N_spacing):N_spacing:round_to_next_N(x_l[2], N=N_spacing)
end