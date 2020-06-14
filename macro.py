"""Compartmental or macro modelling of infectious diseases.

This module provides an API for macro or compartmental models of
infectious disease epidemics.

Users provide model specifications. The output is a time series
showing changes to the models. This output can be converted to
a flat table or exported to a CSV file.

The main functions are:
    * simulate - Takes a list of model specifications and ex
    * simulateSeries - multiple lists of model specifications and
                       runs them in parallel.

Meaningful compartment prefixes
-------------------------------

S - Susceptible  (See delta_S_I and delta_S_I1)
E - Exposed  (See delta_S_I and delta_S_I1)
I - Infectious (See delta_S_I and delta_S_I1)
A - Asymptomatic
T - On treatment (See delta_S_I1 and treatment_infectiousness)
N - Population size. Strictly reserved. Do not prefix any compartment N.
D - Dead (not included in totalling N)
B - Birth (not tallied)
R - Recovered
M - Maternal immunity
V - Vaccinated

"""

import csv
from copy import deepcopy
import random
from multiprocessing import Pool
from typing import List, Dict
import os

Model = Dict
ModelList = List[Model]
ModelListSeries = List[ModelList]


def delta_S_I(from_to, prod, compartments, totals, model=None):
    from_, to_ = from_to.split("_")
    return prod * compartments[from_] * totals[to_] / totals['N']


def delta_X_Y(from_to, prod, compartments, totals, model=None):
    from_, _ = from_to.split("_")
    return prod * compartments[from_]


def delta_birth_X(from_to, prod, compartments, totals, model=None):
    _, to_ = from_to.split("_")
    return prod * compartments[to_]


def make_parameters(dictionary: Dict) -> Dict:
    parameters = deepcopy(PARAMETERS)
    for key, val in dictionary.items():
        if key in parameters:
            if isinstance(parameters[key], dict):
                for key2, val2 in dictionary[key].items():
                    parameters[key][key2] = val2
            else:
                parameters[key] = val
        else:
            parameters[key] = val
    return parameters


def traverse(group: Dict) -> Dict:
    yield group
    if 'groups' in group:
        for group in group['groups']:
            yield from traverse(group)


def sum_infectiousness(model: Model) -> float:
    total = 0.0
    ti = model['parameters']['treatment_infectiousness']
    ai = model['parameters']['asymptomatic_infectiousness']
    for group in traverse(model):
        if 'compartments' in group:
            for key, value in group['compartments'].items():
                if key[0] == 'I':
                    total += value
                if key[0] == 'T':
                    total += ti * value
                if key[0] == 'A':
                    total += ai * value
    return total


def delta_S_I1(from_to, prod, compartments, totals, model=None):
    from_, _ = from_to.split("_")
    infections = 0
    infections = sum_infectiousness(model)
    return prod * compartments[from_] * infections / totals['N']


def reduce_infectivity(model: Model):
    reduction = model['parameters'].get('reduce_infectivity', 1.0)

    for group in traverse(model):
        if 'transitions' in group:
            for key, value in group['transitions'].items():
                from_, to_ = key.split("_")
                if from_[0] == 'S' and (to_[0] == 'E' or to_[0] == 'I'):
                    group['transitions'][key] *= reduction


# These are the default parameters
PARAMETERS = {
    'from': 0,
    'to': 365,
    'record_frequency': 50,  # Write results every X iterations
    # Multiply S_E or S_I by X every iteration if reduce_infectivity
    # function executed
    'reduce_infectivity': 1.0,
    'asymptomatic_infectiousness': 1.0,
    'treatment_infectiousness': 1.0,
    'noise': 0.0,
    'discrete': False,
    'record_first': True,
    'record_last': True,
    'transition_funcs': {
        'S_I': delta_S_I,
        'S_E': delta_S_I,
        'S_I1': delta_S_I1,
        'S_E1': delta_S_I1,
        'B_S': delta_birth_X,
        'default': delta_X_Y
    },
    'before_funcs': [],
    'after_funcs': [],
}


def _update_compartments(model, totals, group=None,
                         transitions=None, parameters=None):
    if group is None:
        group = model

    if transitions is None:
        transitions = {}
    t = transitions.copy()
    if 'transitions' in group:
        t.update(group['transitions'])

    if parameters is None:
        parameters = {}
    if 'parameters' in group:
        parameters.update(group['parameters'])

    if 'compartments' in group:
        compartments = group['compartments']
        deltas = {}
        for key, value in t.items():
            func = None
            if key in parameters['transition_funcs']:
                func = parameters['transition_funcs'][key]
            else:
                func = parameters['transition_funcs']['default']
            val = func(key, value, compartments, totals, model)
            noise = parameters['noise']
            if noise:
                val *= random.uniform(1.0 - noise, 1.0 + noise)
            if parameters['discrete']:
                val = round(val, 0)
            deltas[key] = val
        for key, value in deltas.items():
            from_, to_ = key.split("_")
            compartments[from_] -= value
            compartments[to_] += value

    if 'groups' in group:
        for group in group['groups']:
            _update_compartments(model, totals, group, t, parameters)


def calc_totals(model: Model) -> List[float]:
    totals = {'N': 0}
    for group in traverse(model):
        if 'compartments' in group:
            for key, value in group['compartments'].items():
                if key in totals:
                    totals[key] += value
                else:
                    totals[key] = value
                if key[0] != 'D':
                    totals['N'] += value
    return totals


def _sum_compartment(total_dict, compartment_prefixes):
    total = 0
    for key, value in total_dict.items():
        if key[0] in compartment_prefixes:
            total += value
    return total


def _R0(modelSeries):
    totals = [calc_totals(model) for model in modelSeries]
    infections = [(_sum_compartment(total, {"I", "E"}),
                   _sum_compartment(total, {"S"}),
                   _sum_compartment(total, {"N"})) for total in totals]
    r0_eq_1 = 0
    for i in range(2, len(infections) - 1):
        if infections[i][0] < infections[i-1][0]:
            r0_eq_1 = i - 1
            break
    front = infections[r0_eq_1 - 1][0]
    middle = infections[r0_eq_1][0]
    back = infections[r0_eq_1 + 1][0]
    diff = (middle - front) + (middle - back) / 2.0
    s = infections[r0_eq_1][1] - diff
    r0 = 1 / (s / infections[r0_eq_1][2])
    r1 = 1 / (infections[r0_eq_1][1] / infections[r0_eq_1][2])
    s = infections[r0_eq_1][1] + diff
    r2 = 1 / (s / infections[r0_eq_1][2])
    r3 = 1 / (infections[r0_eq_1 - 1][1] / infections[r0_eq_1][2])
    r4 = 1 / (infections[r0_eq_1 + 1][1] / infections[r0_eq_1][2])
    r5 = 1 / ((infections[r0_eq_1 - 1][1] + diff) / infections[r0_eq_1][2])
    return (r0, r1, r2, r3, r4, r5)


# This needs to be reconceptualised. Highly inaccurate for large R0
def R0(modelListSeries: ModelListSeries) -> List[float]:
    if len(modelListSeries) < 4:
        return None
    r0 = []
    for i in range(len(modelListSeries[0])):
        r0.append(_R0([m[i] for m in modelListSeries]))
    return r0


def _iterate_model(modelList, ident=None):
    modelListSeries = []
    firstModelList = []
    for model in modelList:
        if model['parameters']['record_first']:
            if ident is not None:
                model['ident'] = ident
            model['iteration'] = 0
            firstModelList.append(deepcopy(model))
    if len(firstModelList) > 0:
        modelListSeries.append(firstModelList)
    from_ = min([m['parameters']['from'] for m in modelList])
    to_ = max([m['parameters']['to'] for m in modelList])

    for iteration in range(from_, to_):
        iterationModelList = []
        for model in modelList:
            if iteration < model['parameters']['from'] or \
               iteration >= model['parameters']['to']:
                break
            model['iteration'] = iteration + 1
            if ident is not None:
                model['ident'] = ident
            for func in model['parameters']['before_funcs']:
                func(model)
            totals = calc_totals(model)
            _update_compartments(model, totals)
            for func in model['parameters']['after_funcs']:
                func(model)
            if (iteration + 1) % model['parameters']['record_frequency'] == 0:
                iterationModelList.append(model)
        if len(iterationModelList) > 0:
            modelListSeries.append(deepcopy(iterationModelList))

    lastModelList = []
    for model in modelList:
        if model['parameters']['record_last']:
            if ident is not None:
                model['ident'] = ident
            model['iteration'] = to_
            lastModelList.append(deepcopy(model))
    if len(lastModelList) > 0:
        modelListSeries.append(lastModelList)

    return modelListSeries


def _get_header(model, concat_names=None):
    biggest = []
    current = []
    firstName = True
    i = 0
    for group in traverse(model):
        if 'ident' in group:
            current.append('ident')
        if 'iteration' in group:
            current.append('iter')
        if 'name' in group:
            if concat_names and firstName is False:
                pass
            else:
                if concat_names and firstName:
                    current.append('name')
                else:
                    current.append('name_' + str(i))
                firstName = False
                i += 1
        if 'compartments' in group:
            for key in group['compartments']:
                current.append(key)
            if len(current) > len(biggest):
                biggest = current.copy()
            current = []
    return biggest


def model_to_table(model: Model, concat_names=None) -> List[List]:
    table = []
    current = []
    idents = []
    prevIdents = idents
    names = []
    prevNames = names
    for group in traverse(model):
        if 'ident' in group:
            idents.append(group['ident'])
        if 'iteration' in group:
            idents.append(group['iteration'])
        if 'name' in group:
            names.append(group['name'])
        if 'compartments' in group:
            i = max(len(prevIdents) - len(idents), 0)
            current += prevIdents[:i] + idents

            i = max(len(prevNames) - len(names), 0)
            if concat_names is None:
                current += prevNames[:i] + names
            else:
                current.append(concat_names.join(prevNames[:i] + names))
            for key, val in group['compartments'].items():
                current.append(val)
            table.append(current.copy())
            current = []
            if len(prevIdents) < len(idents):
                prevIdents = idents.copy()
            idents = []
            if len(prevNames) < len(names):
                prevNames = names.copy()
            names = []
    return table


def modelList_to_table(modelList: ModelList, header=True,
                       concat_names=None) -> List[List]:
    table = []
    if header:
        table.append(_get_header(modelList[0], concat_names))
    for model in modelList:
        table += model_to_table(model, concat_names)
    return table


def table_to_csv(table: List[List], csvfile: str, delimiter=',', quotechar='"',
                 quoting=csv.QUOTE_MINIMAL):
    with open(csvfile, 'w', newline='') as csvfile:
        out = csv.writer(csvfile, delimiter=delimiter,
                         quotechar=quotechar, quoting=quoting)
        for row in table:
            out.writerow(row)


def series_to_table(modelListSeries: ModelListSeries, header=True,
                    concat_names=None) -> List[List]:
    table = []
    if header:
        table.append(_get_header(modelListSeries[0][0], concat_names))
    for modelList in modelListSeries:
        tbl = modelList_to_table(modelList, False, concat_names)
        for row in tbl:
            table.append(row)
    return table


def series_to_csv(modelListSeries: ModelListSeries, csvfile: str,
                  header=True, delimiter=',',
                  quotechar='"', quoting=csv.QUOTE_MINIMAL, concat_names=None):
    table = series_to_table(modelListSeries, header, concat_names)
    table_to_csv(table, csvfile, delimiter, quotechar, quoting)


def simulate(modelList: Model, ident=None) -> ModelListSeries:
    results = []
    p = deepcopy(PARAMETERS)
    for model in modelList:
        m = deepcopy(model)
        if 'parameters' in m:
            p.update(m['parameters'])
        m['parameters'] = p
        results.append(m)
    return _iterate_model(results, ident)


def _simulate(m):
    return simulate(m[0], m[1])


def simulate_series(modelListSeries: ModelListSeries,
                    processes=os.cpu_count()) -> ModelListSeries:

    mls_with_ident = [(m[0], m[1]) for m in
                      zip(modelListSeries, range(len(modelListSeries)))]
    with Pool(processes=processes) as pool:
        output = pool.map(_simulate, mls_with_ident)
    results = []
    for r in output:
        results += r
    return results
