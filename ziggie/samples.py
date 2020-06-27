from copy import deepcopy
import random
from ziggie import macro


def mix_models(model, modelList):

    def _calc(val1, val2, prop, noise):
        result = min(val1, val2) * prop
        result *= random.uniform(1.0 - noise, 1.0 + noise)
        return result

    noise = model['parameters']['noise']
    informal = modelList[0]['groups']
    formal = modelList[1]['groups']
    rural = modelList[2]['groups']

    for i in range(len(informal)):
        for key in informal[i]['compartments']:
            if key in ['S', 'E', 'Im', 'A', 'R']:
                delta_f_i = _calc(informal[i]['compartments'][key],
                                  formal[i]['compartments'][key], 0.02,
                                  noise)
                delta_i_r = _calc(informal[i]['compartments'][key],
                                  rural[i]['compartments'][key], 0.01,
                                  noise)
                delta_r_f = _calc(formal[i]['compartments'][key],
                                  rural[i]['compartments'][key], 0.01,
                                  noise)

                if model['iteration'] % 2 == 0:
                    delta_f_i = -delta_f_i
                    delta_i_r = -delta_i_r
                    delta_r_f = -delta_r_f

                formal[i]['compartments'][key] += delta_f_i
                informal[i]['compartments'][key] -= delta_f_i

                informal[i]['compartments'][key] += delta_i_r
                rural[i]['compartments'][key] -= delta_i_r

                rural[i]['compartments'][key] += delta_r_f
                formal[i]['compartments'][key] -= delta_r_f



class MacroModels():

    def simple(self):

        return {
            'name': 'Simple model',
            'compartments': {
                'S': 57000000,
                'I': 1,
                'R': 0
            },
            'transitions': {
                'S_I': 0.6,
                'I_R': 0.1
            },
        }

    def complicated(self):
        return {
            'name': 'Van Wyks Dorp',
            'groups': [
                {
                    'transitions': {
                        'I_R': 0.1,
                    },
                    'groups': [
                        {
                            'name': 'Male',
                            'transitions': {
                                'S_I': 0.3,
                            },
                            'groups': [
                                {
                                    'name': '0-50',
                                    'compartments': {
                                        'S': 290.0,
                                        'I': 1.0,
                                        'R': 0.0,
                                    }
                                },
                                {
                                    'name': '50-100',
                                    'compartments': {
                                        'S': 200.0,
                                        'I': 0.0,
                                        'R': 0.0,
                                    }
                                }
                            ]
                        },
                        {
                            'name': 'Female',
                            'transitions': {
                                'S_I': 0.15,
                            },
                            'groups': [
                                {
                                    'name': '0-50',
                                    'compartments': {
                                        'S': 310.0,
                                        'I': 0.0,
                                        'R': 0.0,
                                    }
                                },
                                {
                                    'name': '50-100',
                                    'compartments': {
                                        'S': 290.0,
                                        'I': 0.0,
                                        'R': 0.0,
                                    }
                                }
                            ]
                        },
                    ]
                }
            ]
        }

    def granich(self, treat_only_I4=False):
        g = {
            'name': 'Granich_et_al_2009',
            'compartments': {
                'S':  30000000.0,
                'I1':   750000.0,
                'I2':   750000.0,
                'I3':   750000.0,
                'I4':   750000.0,
                'T1':   750000.0,
                'T2':   750000.0,
                'T3':   750000.0,
                'T4':   750000.0,
                'D':         0.0,
                'DI':        0.0,
                'DT':        0.0,
                'B':         0.0
            },
            'transitions': {
                'S_I1':  0.001,
                'I1_I2': 1.0/(365.25 * 2.0),
                'I2_I3': 1.0/(365.25 * 2.0),
                'I3_I4': 1.0/(365.25 * 2.0),

                'I1_T1': 0.005,
                'I2_T2': 0.005,
                'I3_T3': 0.01,
                'I4_T4': 0.005,
                'T1_I1': 0.0001,
                'T2_I2': 0.0001,
                'T3_I3': 0.0001,
                'T4_I4': 0.0001,

                'I4_D': 1.0 / (365.25*70),
                'S_D': 1.0 / (365.25*70),
                'I1_D': 1.0 / (365.25*70),
                'I2_D': 1.0 / (365.25*70),
                'I3_D': 1.0 / (365.25*70),
                'T1_D': 1.0 / (365.25*70),
                'T2_D': 1.0 / (365.25*70),
                'T3_D': 1.0 / (365.25*70),
                'T4_D': 1.0 / (365.25*70),
                'I4_DI': 0.005,
                'T4_DT': 0.0001,

                'B_S': (1000000/50000000)/365.25
            },
            'parameters': {
                'to': 365*20,
                'after_funcs': [macro.reduce_infectivity, ],
                'reduce_infectivity': 0.9999,
                'treatment_infectiousness': 0.001,
                'noise': 0.1
            },
        }

        if treat_only_I4:
            g['transitions']['I1_T1'] = 0.0
            g['transitions']['I2_T2'] = 0.0
            g['transitions']['I3_T3'] = 0.0
            g['transitions']['I3_T4'] = 0.01

        return g

    def corona(self):
        parameters = {
            'to': 365,
            'record_frequency': 1,
            'record_last': False,
            'noise': 0.1,
            'asymptomatic_infectiousness': 0.75,
            'reduce_infectivity': 0.999,
            'after_funcs': [macro.reduce_infectivity, ],
            'transition_funcs': {
                'S_E': macro.delta_S_I1,
            }
        }

        parameters_rural = deepcopy(parameters)
        parameters_rural['after_funcs'] = [macro.reduce_infectivity,
                                           mix_models, ]
        return [
            {
                'name': 'Urban informal',
                'parameters': parameters,
                'transitions': {
                    'S_E': 0.4,
                    'E_A': 0.125,
                    'E_Im': 0.125,
                    'Im_Ic': 0.2,
                    'Ic_R': 0.2,
                    'A_R': 0.2
                },
                'groups': [
                    {
                        'name': '0-24',
                        'transitions': {
                            'E_A': 0.25,
                            'E_Im': 0.01,
                            'Ic_D': 0.002
                        },
                        'compartments': {
                            'S': 2100000,
                            'E': 10,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '25-54',
                        'transitions': {
                            'Ic_D': 0.0032
                        },
                        'compartments': {
                            'S': 2100000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '55-',
                        'transitions': {
                            'Ic_D': 0.032
                        },
                        'compartments': {
                            'S': 600000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    }
                ]
            },
            {
                'name': 'Urban formal',
                'parameters': parameters,
                'transitions': {
                    'S_E': 0.3,
                    'E_A': 0.125,
                    'E_Im': 0.125,
                    'Im_Ic': 0.2,
                    'Ic_R': 0.2,
                    'A_R': 0.2
                },
                'groups': [
                    {
                        'name': '0-24',
                        'transitions': {
                            'E_A': 0.25,
                            'E_Im': 0.01,
                            'Ic_D': 0.002
                        },
                        'compartments': {
                            'S': 16940000,
                            'E': 10,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '25-54',
                        'transitions': {
                            'Ic_D': 0.003
                        },
                        'compartments': {
                            'S': 16940000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '55-',
                        'transitions': {
                            'Ic_D': 0.03
                        },
                        'compartments': {
                            'S': 4620000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    }
                ]
            },
            {
                'name': 'Rural',
                'parameters': parameters_rural,
                'transitions': {
                    'S_E': 0.27,
                    'E_A': 0.125,
                    'E_Im': 0.125,
                    'Im_Ic': 0.2,
                    'Ic_R': 0.2,
                    'A_R': 0.2
                },
                'groups': [
                    {
                        'name': '0-24',
                        'transitions': {
                            'E_A': 0.25,
                            'E_Im': 0.01,
                            'Ic_D': 0.002
                        },
                        'compartments': {
                            'S': 7260000,
                            'E': 10,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '25-54',
                        'transitions': {
                            'Ic_D': 0.0035
                        },
                        'compartments': {
                            'S': 7260000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    },
                    {
                        'name': '55-',
                        'transitions': {
                            'Ic_D': 0.035
                        },
                        'compartments': {
                            'S': 2000000,
                            'E': 0,
                            'Im': 0,
                            'Ic': 0,
                            'A': 0,
                            'R': 0,
                            'D': 0
                        }
                    }

                ]
            }
        ]
