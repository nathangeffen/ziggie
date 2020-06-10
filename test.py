import csv
from copy import deepcopy
import macro
import tempfile
import unittest

class TestSimple(unittest.TestCase):

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

    def test_simulate(self):
        res = macro.simulate([self.simple()])
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        self.assertEqual(len(res), 9)
        self.assertEqual(len(res[0]), 1)

        totals = macro.calc_totals(res[-2][0])

        self.assertAlmostEqual(totals['N'], 57000001.0,
                               msg="Sum of compartments")
        table = macro.modelList_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertEqual(table[1], [350, 'Simple model', 65525.67886409864,
                                    7.38650518227673e-07, 56934475.32113517],
                         msg="Table row has correct values")
        (_, filename) = tempfile.mkstemp()
        macro.series_to_csv(res, filename)
        with open(filename, newline='') as csvfile:
            c = csv.reader(csvfile)
            t = [r for r in c]
            self.assertEqual(table[0], t[0])
            for pair in zip(table[-1], t[-2]):
                self.assertEqual(str(pair[0]), str(pair[1]))

    def test_after_func(self):
        parameters = {'before_funcs': [macro.reduce_infectivity, ]}
        parameters['reduce_infectivity'] = 1.0
        simp = self.simple()
        before_S_I = simp['transitions']['S_I']
        simp['parameters'] = parameters
        res = macro.simulate([simp])
        after_S_I = res[-1][0]['transitions']['S_I']
        self.assertAlmostEqual(res[-1][0]['compartments']['S'],
                               65525.67886409458,
                               "Susceptible with reduce_infectivity at 1.0")
        self.assertAlmostEqual(before_S_I, 0.6,
                          "Infectivity reducing function before at 1.0")
        self.assertAlmostEqual(after_S_I, 0.6,
                          "Infectivity reducing function after at 1.0")

        parameters['reduce_infectivity'] = 0.99
        simp = self.simple()
        before_S_I = simp['transitions']['S_I']
        simp['parameters'] = parameters
        res = macro.simulate([simp])
        after_S_I = res[-1][0]['transitions']['S_I']
        self.assertAlmostEqual(res[-1][0]['compartments']['S'],
                               2894902.5260503027,
                               msg="Susceptible with reduce_infectivity at 0.99")
        self.assertAlmostEqual(before_S_I, 0.6,
                               msg="Infectivity reducing function before at 1.0")
        self.assertAlmostEqual(after_S_I, 0.015310778671374667,
                          msg="Infectivity reducing function at 0.99")

class TestComplicated(unittest.TestCase):

    def complicated(self):
        return  {
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

    def test_traverse(self):
        groups = self.complicated()
        generator = macro.traverse(groups)
        lst = list(generator)
        self.assertEqual(len(lst), 8, "Test traverse yields all groups")
        lst2 = [l['name'] for l in lst if 'name' in l]
        self.assertEqual(len(lst2), 7, "Test traverse yields correct groups")

    def test_simulate(self):
        totals = macro.calc_totals(self.complicated())
        self.assertAlmostEqual(totals['N'], 1091.0,
                               msg="Sum of compartments")
        res = macro.simulate([self.complicated()])
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        self.assertEqual(len(res),9)
        self.assertEqual(len(res[0]),1)
        totals = macro.calc_totals(res[-2][0])

        self.assertAlmostEqual(totals['N'], 1091.0,
                               msg="Sum of compartments")
        table = macro.modelList_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'name_1', 'name_2',
                                    'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertEqual(table[1],
                         [350, 'Van Wyks Dorp', 'Male', '0-50',
                          25.59087276881619, 6.371396740780498e-07,
                          265.4091265940443],
                         msg="Table row has correct values")
        (_, filename) = tempfile.mkstemp()
        macro.series_to_csv(res, filename)
        with open(filename, newline='') as csvfile:
            c = csv.reader(csvfile)
            t = [r for r in c]
            self.assertEqual(table[0], t[0])
            for pair in zip(table[-1], t[-5]):
                self.assertEqual(str(pair[0]), str(pair[1]))

        table = macro.modelList_to_table(res[-2], concat_names="|")
        (_, filename) = tempfile.mkstemp()
        macro.series_to_csv(res, filename, concat_names="|")
        self.assertEqual(len(table[0]), 5)
        self.assertEqual(len(table[0]), len(table[-1]))
        with open(filename, newline='') as csvfile:
            c = csv.reader(csvfile)
            t = [r for r in c]
            self.assertEqual(table[0], t[0])
            for pair in zip(table[-1], t[-5]):
                self.assertEqual(str(pair[0]), str(pair[1]))

class TestNoise(unittest.TestCase):

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
            'parameters': {
                'noise': 0.05
            }
        }

    def test_simulate(self):
        res = macro.simulate([self.simple()])
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        totals = macro.calc_totals(res[-1][0])

        self.assertAlmostEqual(totals['N'], 57000001.0, 1,
                               msg="Sum of compartments")
        table = macro.modelList_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertGreater(table[1][2], 50000.0,
                           msg="Random fluctuations give reasonable result")
        self.assertGreater(80000.0, table[1][2],
                           msg="Random fluctuations give reasonable result")

class TestGranich(unittest.TestCase):

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
            'parameters': macro.make_parameters({
                'to': 365*20,
                'after_funcs':[macro.reduce_infectivity,],
                'reduce_infectivity': 0.9999,
                'treatment_infectiousness': 0.001,
                'noise': 0.1
            }),
        }

        if treat_only_I4:
            g['transitions']['I1_T1'] = 0.0
            g['transitions']['I2_T2'] = 0.0
            g['transitions']['I3_T3'] = 0.0
            g['transitions']['I3_T4'] = 0.01

        return g

    def test_simulate(self):
        treat_all = self.granich(False)
        treat_deferred = self.granich(True)
        results_treat = macro.simulate([treat_all])
        results_deferred = macro.simulate([treat_deferred])

        totals_before_t = macro.calc_totals(results_treat[0][0])
        totals_after_t = macro.calc_totals(results_treat[-1][0])
        n_before = sum([value for key,value in totals_before_t.items()
                        if key != 'N'])
        n_after = sum([value for key, value in totals_after_t.items()
                       if key != 'N'])
        self.assertAlmostEqual(n_before, n_after, 6)
        self.assertAlmostEqual(n_after, 36000000, 6)

        totals_before_d = macro.calc_totals(results_deferred[0][0])
        totals_after_d = macro.calc_totals(results_deferred[-1][0])
        n_before = sum([value for key,value in totals_before_d.items()
                        if key != 'N'])
        n_after = sum([value for key, value in totals_after_d.items()
                       if key != 'N'])
        self.assertAlmostEqual(n_before, n_after, 6)
        self.assertAlmostEqual(n_after, 36000000, 6)

        self.assertGreater(totals_after_t['DI'], 1000)
        self.assertGreater(totals_after_t['DT'], 1000)
        self.assertGreater(totals_after_d['DI'], totals_after_t['DI'])
        self.assertGreater(totals_after_d['DT'], totals_after_t['DT'])
        self.assertGreater(totals_after_t['N'], totals_after_d['N'])
        self.assertGreater(-totals_after_t['B'], -totals_after_d['B'])

if __name__ == '__main__':
    unittest.main()
