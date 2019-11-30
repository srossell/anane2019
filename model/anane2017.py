"""
Kinetic model of overflow metabolism in Eschericia coli by acetate cycling
Anane, Lopez, Neubauer & Cruz Bournazou (2017)
Biochemical Engineering Journal 125: 23-30
"""

id_sp = [
        'V',  #       = y(1); volume
        'X',  #       = y(2); biomass
        'S',  #       = y(3); substrate
        'A',  #       = y(4); acetate
        'DOT',  #     = y(5); dissolved oxygen, no probe response
        'DOTm',  #    = y(6); dissolved oxygen, with probe response time
        # 'F',  #       = y(7); feed
        ]

params = {
        'Kap'  : 0.5052, # 0.5088
        'Ksa'  : 0.0134, # 0.0128
        'Ko'   : 0.0001, # 0.0001
        'Ks'   : 0.0370, # 0.0381
        'Kia'  : 1.2399, # 1.2602
        'Kis'  : 2.1231, # 1.8383
        'pAmax': 0.2268, # 0.2286
        'qAmax': 0.1148, # 0.1148
        'qm'   : 0.0129, # 0.0133
        'qSmax': 0.6356, # 0.6350
        'Yas'  : 0.9097, # 0.8938
        'Yoa'  : 0.5440, # 0.5221
        'Yxa'  : 0.5718, # 0.5794
        'Yem'  : 0.5333, # 0.5321
        'Yos'  : 1.5620, # 1.5722
        'Yxsof': 0.2268, # 0.2290
        # Fermentation constants
        'Cs': 0.391,  # carbon content of glucose (ref:enfors)
        'Cx': 0.488,  # carbon content of biomas(e coli)
        'H' : 14000,  # Henry's constant to convert DO from mol to % DO
        # INPUTS
        'Si': 300,
        'mufeed': 0.222,
        'DOTstar': 99,
        'Kla': 220,
        'tau': 35,
        'F': 0,
        }

functions = {
        'qS':  '(qSmax/(1+(A/Kia)))*(S/(S+Ks))',
        'qSof':'pAmax*(qS/(qS+Kap))',
        'pA':  'qSof*Yas',
        'qSox':'(qS-qSof)*(DOT/(DOT+Ko))',
        'qSan':'(qSox-qm)*Yem*Cx/Cs',
        'qsA': '(qAmax/(1+(qS/Kis)))*(A/(A+Ksa))',
        'qA':  ' pA-qsA',
        'mu':  '(qSox-qm)*Yem + qsA*Yxa + qSof*Yxsof',
        'qO':  'Yos*(qSox-qSan)+qsA*Yoa',
        'Kp':  '(1/tau)*3600',
        }

id_rs = ['feed', 'v1', 'v2', 'v3', 'v4', 'v5']

rates = {
        'feed': 'F',
        'v1': 'X*(mu - (F/V))',
        'v2': 'F/V*(Si - S) - (qS*X)',
        'v3': 'qA*X-(F/V)*A',
        'v4': 'Kla*(DOTstar-DOT)-qO*X*H',
        'v5': 'Kp*(DOT-DOTm)',
        }

mass_balances = {
        'V': {'feed': 1},
        'X': {'v1': 1},
        'S': {'v2': 1},
        'A': {'v3': 1},
        'DOT': {'v4': 1},
        'DOTm': {'v5': 1},
        }

y0 = [2, 0.17, 4.94, 0.0129, 98, 98]

if __name__ == '__main__':
    from scipy.integrate import solve_ivp
    import matplotlib.pyplot as plt
    import numpy as np

    from bcodes.ratevector import create_rate_vector
    from bcodes.ratevector import subs_id_by_value
    from bcodes.stoichiometrymatrix import build_stoichiometry_matrix

    ####################
    # From function strings to lambda functions in the global scope
    func_trans_dict = dict(
            **params,
            **dict(zip(
                id_sp,
                ['y[{}]'.format(i) for i, sp in enumerate(id_sp)]
                )),
            **{
                'qS': 'qS(y)',
                'qSof': 'qSof(y)',
                'qSox': 'qSox(y)',
                'qsA': 'qsA(y)',
                'pA': 'pA(y)',
                'qSan': 'qSan(y)',
                }
            )

    # Functions to execute
    functions2exec = {}
    for f in functions:
        functions2exec[f] = subs_id_by_value(functions[f], func_trans_dict)
    # write functions
    for f in functions2exec:
        exec('{} = lambda y: {}'.format(f, functions2exec[f]))

    # Building ODEs
    v = create_rate_vector(id_sp, id_rs, rates,
            dict(
                **params,
                **{
                     'qS': 'qS(x)',
                     'qA': 'qA(x)',
                     'qO': 'qO(x)',
                     'Kp': 'Kp(x)',
                     'mu': 'mu(x)'
                     }
                ),
            eval_scope = {
                'qS': qS,
                'qSof': qSof,
                'pA': pA,
                'qSox': qSox,
                'qSan': qSan,
                'qA': qA,
                'mu': mu,
                'qO': qO,
                'Kp': Kp
                }
            )

    S = build_stoichiometry_matrix(id_sp, id_rs, mass_balances)

    def odes(t, y):
        return np.dot(S, v(y))

    sol = solve_ivp(odes, [0, 10], y0)

    # Plotting
    fig, ax = plt.subplots(ncols=2, nrows=3)
    for i, sp in enumerate(id_sp):
        ax.flatten()[i].plot(sol.t, sol.y.T[:, i])
        ax.flatten()[i].set_title(sp)
    fig.tight_layout()
    plt.show()




















