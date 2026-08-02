using Test
using RichardsonExtrapolationUQ


r = [45.3028678412486, 46.306074246781954, 46.74593039574231, 46.865808943526055]
h = [0.0078125, 0.00390625, 0.001953125, 0.0009765625]

e = RichardsonExtrapolationUQ.richextrapol_uq(r, h; W = 2/3)
@test e.q_star == 46.923362280076304
@test e.q_star_ci == 0.09041059858140653
@test e.beta_star == 1.735905743595649
@test e.beta_star_ci == 0.3980337138577599
@test e.success

@info "q_ci = $(e.q_star) +/- $(e.q_star_ci), beta_ci = $(e.beta_star) +/- $(e.beta_star_ci)"


r = [45.5, 46.1, 46.81, 46.8]
h = [0.0078125, 0.00390625, 0.001953125, 0.0009765625]

e = RichardsonExtrapolationUQ.richextrapol_uq(r, h; W = 2/3)
@test e.success == false


