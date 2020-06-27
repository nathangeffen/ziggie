"""Tests of macro.py module.

This module runs a series of tests on the macro API. It should be executed
every time code is committed.

It can be executed thus:

python test.py

Requires Python 3.5 or higher.
"""

import csv
from ziggie import macro, samples
import tempfile
import unittest


class TestSimple(unittest.TestCase):

    def simple(self):

        return samples.MacroModels().simple()

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
                               msg="Susceptible with reduce_infectivity==0.99")
        self.assertAlmostEqual(before_S_I, 0.6,
                               msg="Infectivity reduction before")
        self.assertAlmostEqual(after_S_I, 0.015310778671374667,
                               msg="Infectivity reduction after")


class TestComplicated(unittest.TestCase):

    def complicated(self):
        return samples.MacroModels().complicated()

    def test_traverse(self):
        groups = self.complicated()
        generator = macro.traverse(groups)
        lst = list(generator)
        self.assertEqual(len(lst), 8, "Test traverse yields all groups")
        lst2 = [l1['name'] for l1 in lst if 'name' in l1]
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
        self.assertEqual(len(res), 9)
        self.assertEqual(len(res[0]), 1)
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
        return samples.MacroModels().granich(treat_only_I4)

    def test_simulate(self):
        treat_all = self.granich(False)
        treat_deferred = self.granich(True)
        results_treat = macro.simulate([treat_all])
        results_deferred = macro.simulate([treat_deferred])

        totals_before_t = macro.calc_totals(results_treat[0][0])
        totals_after_t = macro.calc_totals(results_treat[-1][0])
        n_before = sum([value for key, value in totals_before_t.items()
                        if key != 'N'])
        n_after = sum([value for key, value in totals_after_t.items()
                       if key != 'N'])
        self.assertAlmostEqual(n_before, n_after, 6)
        self.assertAlmostEqual(n_after, 36000000, 6)

        totals_before_d = macro.calc_totals(results_deferred[0][0])
        totals_after_d = macro.calc_totals(results_deferred[-1][0])
        n_before = sum([value for key, value in totals_before_d.items()
                        if key != 'N'])
        n_after = sum([value for key, value in totals_after_d.items()
                       if key != 'N'])
        self.assertAlmostEqual(n_before, n_after, 5)
        self.assertAlmostEqual(n_after, 36000000, 5)

        self.assertGreater(totals_after_t['DI'], 1000)
        self.assertGreater(totals_after_t['DT'], 1000)
        self.assertGreater(totals_after_d['DI'], totals_after_t['DI'])
        self.assertGreater(totals_after_d['DT'], totals_after_t['DT'])
        self.assertGreater(totals_after_t['N'], totals_after_d['N'])
        self.assertGreater(-totals_after_t['B'], -totals_after_d['B'])


class TestCorona(unittest.TestCase):

    def corona(self):
        return samples.MacroModels().corona()

    def check_results(self, results):
        self.assertEqual(len(results), 366)
        Nb = macro.grand_sum_totals([macro.calc_totals(m)
                                     for m in results[0]])
        Na = macro.grand_sum_totals([macro.calc_totals(m)
                                     for m in results[-1]])
        self.assertAlmostEqual(Nb, Na, 5)
        totals = macro.sum_totals([macro.calc_totals(m)
                                   for m in results[-1]])
        self.assertGreater(totals['D'],  100000)
        self.assertLess(totals['D'], 160000)
        # Check that informal settlement > urban > rural
        proportions = [
            t['D'] / t['N'] for t in
            [macro.calc_totals(m) for m in results[-1]]
        ]

        self.assertTrue(proportions[0] > proportions[1])
        self.assertTrue(proportions[1] > proportions[2])
        self.assertGreater(proportions[2], 0.00015)

    def test_simulate(self):
        results = macro.simulate(self.corona())
        self.assertEqual(len(results), 366)
        self.check_results(results)

    def test_parallel(self):
        resultSeries = macro.simulate_series([TestCorona().corona()
                                              for _ in range(10)])
        self.assertEqual(len(resultSeries), 3660)
        self.assertEqual(len(resultSeries[0]), 3)
        self.assertEqual(len(resultSeries[3659]), 3)
        results = [r for r in resultSeries if r[0]['iteration'] == 365]

        N_0 = macro.calc_totals(results[0][0])['N']
        N_1 = macro.calc_totals(results[1][0])['N']
        self.assertNotEqual(N_0, N_1)
        self.assertLess(abs(N_0 - N_1), 180000)

        for i in range(10):
            results = [r for r in resultSeries if r[0]['ident'] == i]
            self.check_results(results)
            results = [r for r in resultSeries if r[1]['ident'] == i]
            self.check_results(results)


if __name__ == '__main__':
    unittest.main()
