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

    def test_model(self):
        res = macro.model(self.simple())
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        totals = macro.calc_totals(res[-2][0])

        self.assertAlmostEqual(totals['N'], 57000001.0,
                               msg="Sum of compartments")
        table = macro.groups_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertEqual(table[1], [350, 'Simple model', 65525.67886409864,
                                    7.38650518227673e-07, 56934475.32113517],
                         msg="Table row has correct values")
        (_, filename) = tempfile.mkstemp()
        macro.results_to_csv(res, filename)
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
        res = macro.model(simp, parameters)
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
        res = macro.model(simp, parameters)
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

    def test_model(self):
        totals = macro.calc_totals(self.complicated())
        self.assertAlmostEqual(totals['N'], 1091.0,
                               msg="Sum of compartments")
        res = macro.model(self.complicated())
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        totals = macro.calc_totals(res[-2][0])

        self.assertAlmostEqual(totals['N'], 1091.0,
                               msg="Sum of compartments")
        table = macro.groups_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'name_1', 'name_2',
                                    'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertEqual(table[1],
                         [350, 'Van Wyks Dorp', 'Male', '0-50',
                          25.59087276881619, 6.371396740780498e-07,
                          265.4091265940443],
                         msg="Table row has correct values")
        (_, filename) = tempfile.mkstemp()
        macro.results_to_csv(res, filename)
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

    def test_model(self):
        res = macro.model(self.simple())
        self.assertEqual(res[0][0]['iteration'], 0)
        self.assertEqual(res[1][0]['iteration'], 50)
        self.assertEqual(res[-1][0]['iteration'], 365)
        self.assertEqual(res[-2][0]['iteration'], 350)
        totals = macro.calc_totals(res[-1][0])

        self.assertAlmostEqual(totals['N'], 57000001.0, 1,
                               msg="Sum of compartments")
        table = macro.groups_to_table(res[-2], True)
        self.assertEqual(table[0], ['iter', 'name_0', 'S', 'I', 'R'],
                         msg="Table header has correct column names")
        self.assertGreater(table[1][2], 50000.0,
                           msg="Random fluctuations give reasonable result")
        self.assertGreater(80000.0, table[1][2],
                           msg="Random fluctuations give reasonable result")

class TestGranich(unittest.TestCase):

    def granich(self):
        return {

        }
if __name__ == '__main__':
    unittest.main()
