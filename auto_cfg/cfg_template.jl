template = Dict([
    ("seed", 0),
    ("d_fitness", 1),
    ("n_population", 100),
    ("n_elite", 10),
    ("n_gen", 50000),
    ("log_gen", 100),
    ("save_gen", 1000),
    ("m_rate", 0.5),
    ("rows", 1),
    ("columns", 25),
    ("recur", 0.9),
    ("n_out", 1),
    ("n_in", 12),
    ("out_m_rate", 0.9),
    ("functions", [
        "f_subtract",
        "f_add",
        "f_div",
        "f_mult",
        "f_sqrt",
        "f_abs",
        "f_exp",
        "f_pow",
        "f_sin",
        "f_cos", 
        "f_tanh", 
        "f_sqrt_xy",
        "f_lt",
        "f_gt"
        ]
    )
])
