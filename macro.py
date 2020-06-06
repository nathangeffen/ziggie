import csv
from copy import deepcopy
import random
from multiprocessing import Pool
import os

def delta_S_I(from_to, prod, compartments, totals):
    return prod * compartments['S'] * totals['I'] / totals['N']

def delta_X_Y(from_to, prod, compartments, totals):
    return prod * compartments[from_to[0]]

def delta_birth_X(from_to, prod, compartments, totals):
    return prod * compartments[from_to[2]]

PARAMETERS = {
    'from': 0,
    'to': 365,
    'record_frequency': 50,
    'reduce_infectivity': 0.999,
    'noise': 0.0,
    'discrete': False,
    'record_first': True,
    'record_last': True,
    'transition_funcs' : {
        'S_I': delta_S_I,
        'S_E': delta_S_I,
        'E_I': delta_X_Y,
        'I_R': delta_X_Y,
        'S_M': delta_X_Y,
        'S_D': delta_X_Y,
        'I_M': delta_X_Y,
        'I_D': delta_X_Y,
        'E_M': delta_X_Y,
        'E_D': delta_X_Y,
        'R_M': delta_X_Y,
        'R_D': delta_X_Y,
        'M_D': delta_X_Y,
        'M_S': delta_X_Y,
        'R_S': delta_X_Y,
        'B_S': delta_birth_X,
        'default': delta_X_Y
    },
    'before_funcs': [],
    'after_funcs': [],
}

TOTALS_TEMPLATE = {
    'S': 0.0, # Susceptible
    'E': 0.0, # Exposed
    'I': 0.0, # Infectious
    'R': 0.0, # Recovered
    'M': 0.0, # Maternal immunity
    'V': 0.0, # Vaccinated
    'D': 0.0, # Dead
    'N': 0.0, # total of all the above less dead
    # 'B': 0.0 Births - but these are not explicitly tracked
}


def traverse(groups):
    if type(groups) == type({}):
        _groups = [groups]
    else:
        _groups = groups

    for group in _groups:
        yield group
        if 'groups' in group:
            yield from traverse(group['groups'])

def reduce_infectivity(groups, parameters):
    for group in traverse(groups):
        if 'transitions' in group:
            for key, value in group['transitions'].items():
                from_, to_ = key.split("_")
                if from_[0] == 'S' and (to_[0] == 'E' or to_[0] == 'I'):
                    group['transitions'][key] *= parameters['reduce_infectivity']

def _update_compartments(group, parameters, totals, transitions):
    inherited_transitions = deepcopy(transitions)
    if 'transitions' in group:
        inherited_transitions.update(group['transitions'])

    if 'compartments' in group:
        compartments = group['compartments']
        deltas = {}
        for key, value in inherited_transitions.items():
            func = None
            if key in parameters['transition_funcs']:
                func = parameters['transition_funcs'][key]
            else:
                func = parameters['transition_funcs']['default']
            val = func(key, inherited_transitions[key], compartments, totals)

            noise = parameters['noise']
            if noise:
                val *= random.uniform(1.0 - noise, 1.0 + noise)
            if parameters['discrete']:
                val = round(val, 0)
            deltas[key] = val
        for key, value in deltas.items():
            from_, to_ = key.split("_")
            if from_[0] is not 'B':
                compartments[from_] -= value
            compartments[to_] += value
    elif 'groups' in group:
        for subgroup in group['groups']:
            _update_compartments(subgroup, parameters, totals,
                                inherited_transitions)

def _calc_group_totals(group, totals):
    if 'compartments' in group:
        compartments = group['compartments']
        for key in totals:
            if key is not 'N' and key in compartments:
                totals[key] += compartments[key]
                if key is not 'D':
                    totals['N'] += compartments[key]
    elif 'groups' in group:
        for subgroup in group['groups']:
            _calc_group_totals(subgroup, totals)

def calc_totals(group, parameters=None, name=None, totals=None,
                transitions=None):
    if totals is None:
        totals = deepcopy(TOTALS_TEMPLATE)
    if parameters is None:
        parameters = deepcopy(PARAMETERS)

    if 'transitions' in group:
        if transitions is None:
            transitions = group['transitions']
        else:
            transitions.update(group['transitions'])

    if name and group.get('name') != name:
        if "groups" in group:
            for subgroup in group["groups"]:
                totals = calc_totals(subgroup, parameters, name, totals,
                                     transitions)
    else:
        _calc_group_totals(group, totals)
    return totals

def _sum_compartment(total_dict, compartment_prefixes):
    total = 0
    for key, value in total_dict.items():
        if key[0] in compartment_prefixes:
            total += value
    return total

def R0(groups):
    if len(groups) < 4:
        return None
    totals = [calc_totals(group) for group in groups]
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
    r5 = 1 / ( (infections[r0_eq_1 - 1][1] + diff) / infections[r0_eq_1][2])
    return [r0, r1, r2, r3, r4, r5]

def _iterate_model(groups, parameters, ident=None):
    arr = []
    if parameters['record_first']:
        for group in groups:
            if ident is not None:
                group['ident'] = ident
            group['iteration'] = 0
        arr.append(deepcopy(groups))
    for iteration in range(parameters['from'], parameters['to']):
        for func in parameters['before_funcs']:
            func(groups, parameters)
        for group in groups:
            totals = calc_totals(group, parameters)
            _update_compartments(group, parameters, totals, {})
            group['iteration'] = iteration + 1
            if ident is not None:
                group['ident'] = ident
        for func in parameters['after_funcs']:
            func(groups, parameters)
        if (iteration + 1) % parameters['record_frequency'] == 0:
            arr.append(deepcopy(groups))
    if parameters['record_last'] and \
       (parameters['to'] % parameters['record_frequency'] != 0):
        arr.append(deepcopy(groups))
    return arr

def _group_to_header(group, table, current, depth=0, run=None):
    if 'ident' in group:
        current.append('ident')
    if run and depth == 0:
        current.append("run")
    if 'iteration' in group:
        current.append('iter')
    if 'name' in group:
        current.append('name_' + str(depth))
        depth += 1
    if 'groups' in group:
        return _group_to_header(group['groups'][0], table, current, depth)
    if 'compartments' in group:
        for key in group['compartments']:
            current.append(key)
        return current


def _group_to_table(group, table, current, run=None):
    if 'ident' in group:
        current.append(group['ident'])
    if run:
        current.append(run)
    if 'iteration' in group:
        current.append(group['iteration'])
    if 'name' in group:
        current.append(group['name'])

    if 'compartments' in group:
        for key, val in group['compartments'].items():
            current.append(val)
        table.append(current.copy())
        current = []

    if 'groups' in group:
        for g in group['groups']:
            _group_to_table(g, table, current.copy(), None)

def groups_to_table(groups, header=True, run=None):
    table = []
    if header:
        table.append(_group_to_header(groups[0], table, [], 0, run))
    for group in groups:
        _group_to_table(group, table, [], run)
    return table

def table_to_csv(table, csvfile, delimiter=',', quotechar='"',
                 quoting=csv.QUOTE_MINIMAL):
    with open(csvfile, 'w', newline='') as csvfile:
        out = csv.writer(csvfile, delimiter=delimiter,
                         quotechar=quotechar, quoting=quoting)
        for row in table:
            out.writerow(row)

def groups_to_csv(groups, csvfile, header=True, delimiter=',', quotechar='"',
                  quoting=csv.QUOTE_MINIMAL, run=None):
    table = groups_to_table(groups, header, run)
    table_to_csv(table, csvfile, delimiter, quotechar, quoting)

def results_to_table(results, header=True, run=None):
    table = []
    if header:
        table.append(_group_to_header(results[0][0], table, [], 0, run))
    for result in results:
        tbl = groups_to_table(result, False, run)
        for row in tbl:
            table.append(row)
    return table

def results_to_csv(results, csvfile, header=True, delimiter=',',
                   quotechar='"', quoting=csv.QUOTE_MINIMAL, run=None):
    if type(results[0][0]) == type([]): # check for result of parallel_model
        result2 = []
        for r in results:
            result2 += r
        results = result2
    table = results_to_table(results, header, run)
    table_to_csv(table, csvfile, delimiter, quotechar, quoting)

def model(groups, parameters=None, ident=None):
    results = []
    if type(groups) == type([]):
        for group in groups:
            results.append(deepcopy(groups))
    elif type(groups) == type({}):
        results.append(deepcopy(groups))
    else:
        raise TypeError("First parameter must be dict or array")

    p = PARAMETERS.copy()
    if parameters:
        p.update(parameters)
    if 'parameters' in results[0]:
        p.update(results[0]['parameters'])

    return _iterate_model(results, p, ident)

def _parallel_model(m):

    return model(m[0], m[1], m[2])

def parallel_model(model):

    parms = [(groups, parameters, i) for (groups, parameters), i in
             zip(models, range(len(models)))]

    with Pool() as p:
        results = p.map(_parallel_model, parms)

    return results
