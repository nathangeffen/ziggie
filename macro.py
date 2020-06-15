"""Compartmental or macro modelling of infectious diseases.

This module provides an API for macro or compartmental models of
infectious disease epidemics.

Users provide model specifications. The output is a time series
showing changes to the models. This output can be converted to
a flat table or exported to a CSV file.

The main functions are:
    * simulate - Takes a list of model specifications and executes it.
    * simulate_series - Multiple lists of model specifications and
                        runs them in parallel.
    * series_to_table - Takes the output of the above and
                        converts to flat list of lists, with
                        each entry representing output for one iteration
                        of a group.
    * series_to_csv - Same as above except outputs a csv file.

A model specification consists of groups. A model is itself the
top-level group.

Each group can itself contain groups.

A group specification is a dictionary containing the following
one or more of the following fields:

    * name - Name of group.
    * transitions - Dictionary of probabilities of
                    transitioning from one compartment to another.
    * compartments - Dictionary of the number of people per
                     compartment. A group must either have a
                     "compartments" entry or "groups" entry but
                     not both.
    * groups: - Array of groups within this group.
    * parameters: - Hyper-parameters for the model. These are usually
                    only specified in the top-level (model) group. If
                    not specified then the default PARAMETERS
                    are copied and used in the model.
                    If the hyper-parameters are partially specified,
                    the unspecified parameters are copied from PARAMETERS.

Here are two examples:

This is a very simple model.

{
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

A ModelList is a list of models. The "simulate" function takes
a list of models as its first parameter. Why? Let's say you have
distinct regions you'd like to model, say London and Manchester.
You would specify a model for each and then put them in a list (i.e. you
now have a ModelList data structure). You can then run "simulate"
on this ModelList.

The output is two distinct sets of independent new model specification
with updated compartments on each iteration.  i.e. a list of ModelLists.

This is called a ModelListSeries.

Here's an example of a ModelList of Covid-19 for South Africa. It is
uncalibrated and for illustrative purposes only.

The ModelList consists of three models: one for urban informal settlements, one
for urban formal suburbs and one for rural areas.  Each of these has groups for
people aged 0-24, 25-54 and 55 or older.

[
    {
        'name': 'Urban informal',
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


If the output of each of the models is independent, what's the use of putting
them in a ModelList? Well, you could specify a custom function in the
"after_funcs" hyper parameter of one of the models. This function could migrate
individuals between the two models or connect them in some other way. See the
TestCorona class in test.py for an example.


Transitions
-----------

A transition name consists of two compartment names separated by
an underscore. The default PARAMETERS consist of a dictionary called
'transition_funcs' which you can override in the model's parameters.
The transition function to execute is looked up in this dictionary. If
it's not found the function specified by 'default' is executed.

For most transitions, a proportion of the "from" compartment is moved
to the "to" compartment. The delta_X_Y function takes care of this.

But for new infections you almost always want
a standard SIR-like equations such as this:

delta = susceptibles * number of contacts per iteration *
                       risk of infection per contact *
                       total number of infectious individuals /
                       total population

susceptibles -= delta
infectious (or exposed) individuals += delta

Two functions are provided to deal with this: the very simple delta_S_I
and the more sophisticated (but slower) delta_S_I1.

Some of the compartment name prefixes are meaningful. A compartment
name generally starts with one of these meaningful prefixes and then
a unique identifier. E.g. I1, I2, I3 or I4 for various stages of infection.

Meaningful compartment prefixes
-------------------------------

S - Susceptible  (See delta_S_I and delta_S_I1)
E - Exposed  (See delta_S_I and delta_S_I1)
I - Infectious (See delta_S_I and delta_S_I1)
A - Asymptomatic (See delta_S_I1 and asymptomatic_infectiousness)
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
from typing import List, Dict, Generator
import os

Model = Dict
Group = Dict
ModelList = List[Model]
ModelListSeries = List[ModelList]


def delta_S_I(from_to, beta, compartments, totals, model=None):
    """Return number of new infections.

    Parameters:
    from_to (str): transition name consisting of two compartment names
                    separated by an underscore (e.g. S_Ic)
    beta (float): effective contact rate
    compartments (dict): dictionary of compartments including the two
                         specified in from_to
    totals (dict): dictionary containing a key 'N' that is the sum of the
                   total population for this model.
    model (Model): Unused but part of function signature
    """
    from_, to_ = from_to.split("_")
    return beta * compartments[from_] * totals[to_] / totals['N']


def delta_X_Y(from_to, prop, compartments, totals, model=None):
    """Return number individuals to be moved from one compartment to another.

    Parameters:
    from_to (str): transition name consisting of two compartment names
                    separated by an underscore (e.g. I_R)
    prop (float): proportion of "from" compartment to move
    compartments (dict): dictionary of compartments including the two
                         specified in from_to
    totals (dict): Unused but part of function signature
    model (Model): Unused but part of function signature
    """
    from_, _ = from_to.split("_")
    return prop * compartments[from_]


def delta_birth_X(from_to, prop, compartments, totals, model=None):
    """Return number new individuals in a population.

    Parameters:
    from_to (str): transition name consisting of a compartment name
                   starting with a B followed by an underscore
                   then another compartment, typically B_S (for
                   number of births in susceptible population)
    prop (float): proportion of "to" compartment that will be added
    compartments (dict): dictionary of compartments including the "from"
                         part of the from_to
    totals (dict): Unused but part of function signature
    model (Model): Unused but part of function signature
    """
    _, to_ = from_to.split("_")
    return prop * compartments[to_]


def _make_parameters(dictionary: Dict) -> Dict:
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


def traverse(group: Group) -> Generator[Group, None, None]:
    """Generate all the groups of a model or group.

    Parameters:
    group: the group whose subgroups to traverse
    """
    yield group
    if 'groups' in group:
        for group in group['groups']:
            yield from traverse(group)


def sum_infectiousness(model: Model) -> float:
    """Return weighted sum of all infection compartments.

    The infection compartments are all those that begin with an 'I'
    (infectious), 'A' (asymptomatic) or 'T' (treatment).
    Each 'A' compartment is multipled by the parameter
    asymptomatic_infectiousness. Each 'T' compartment is multiplied by
    the parameter treatment_infectiousness.

    Parameters:
    model: the model to sum
    """
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


def delta_S_I1(from_to, beta, compartments, totals, model=None):
    """Return number of new infections.

    In contrast to delta_S_I this function calculates the weighted
    infectiousness of the population. In other words it calculates::

    number susceptible * effective contact rate per iteration *
    sum_infectiousness(model) / total population

    Parameters:
    from_to (str): transition name consisting of two compartment names
                    separated by an underscore (e.g. S_I1)
    beta (float): effective contact rate per iteration
    compartments (dict): dictionary of compartments including the two
                         specified in from_to
    totals (dict): dictionary containing a key 'N' that is the sum of the
                   total population for this model.
    model (Model): model to calculate total infectiousness for
    """
    from_, _ = from_to.split("_")
    infections = 0
    infections = sum_infectiousness(model)
    return beta * compartments[from_] * infections / totals['N']


def reduce_infectivity(model: Model, modelList=None):
    """Reduce the values of the effective contact rates of a model.

    This is an optional function to be executed before or after each
    iteration. To use it, add it to the after_funcs or before_funcs
    parameters. It traverses the model aand reduces the value of
    all parameters of the S(.*)_[(E|I)(.*) by multiplying them
    by the 'reduce_infectivity' parameter.
    This is meant to simulate the fact that in heterogeneous
    populations, the most susceptible will be infected first.

    Parameters
    model (Model): the model to apply the reductions to
    modelList (ModelList): Unused but part of function signature

    """
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


def calc_totals(model: Model) -> Dict[str, float]:
    """Calculate sum of each compartment in all groups and return dict.

    This function traverses a model and calculates the sum of each compartment
    across all the groups. It also calculates N, the total living population.
    """
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


def grand_sum_totals(totals: List[Dict[str, float]],
                     ignore=['B', 'N', ]) -> float:
    """Calculate the sum of all compartments.

    Takes the output list generated using calc_totals and calculates the sum of
    all the compartments.

    E.g. print(grand_sum_totals([calc_totals[m] for m in model_list])

    Parameters:
    totals (output of calc_totals): list of compartment totals across models
    ignore (list of str): list of compartments to ignore

    """
    total = 0.0
    for t in totals:
        for key, value in t.items():
            if key not in ignore:
                total += value
    return total


def sum_totals(totals: List[Dict[str, float]]) -> List[Dict[str, float]]:
    """Calculate the sum of compartments across multiple models.

    Takes the output list generated using calc_totals and calculates the sum of
    each the compartments.

    E.g. print(sum_totals([calc_totals[m] for m in model_list])

    Parameters:
    totals (output of calc_totals): list of compartment totals across models
    """
    result = {}
    for t in totals:
        for key, value in t.items():
            result[key] = result.get(key, 0) + value
    return result


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


def R0(modelListSeries: ModelListSeries) -> List[float]:
    """Calculate R0 for a time series of model outputs.

    Highly inaccurate and needs to be reconceptualised. Don't use for now.
    """
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
                func(model, modelList)
            totals = calc_totals(model)
            _update_compartments(model, totals)
            for func in model['parameters']['after_funcs']:
                func(model, modelList)
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
    """Create and return a table from a model.

    This converts a model to list of lists, where each list corresponds
    to a group in the model. This is also an interim step to converting
    a model output to a csv file.

    Parameters
    model (Model): model to convert
    concat_names (str): if not None then group names are concatenated,
                        separated by this string
    """
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
    """Create and return a table from a list of models.

    This function "flattens" a list of models so that each row in the
    table corresponds to a group in each of the models in the list.
    It's especially useful for analysing the model output, for
    example with numpy. It is also an interim step to converting a
    list of models to a comma-separated file.

    Parameters
    modelList (ModelList): list of models to convert
    header (bool): whether or not the first row should be a header
    concat_names (str): if not None then group names are concatenated,
                        separated by this string
    """
    table = []
    if header:
        table.append(_get_header(modelList[0], concat_names))
    for model in modelList:
        table += model_to_table(model, concat_names)
    return table


def table_to_csv(table: List[List], csvfile: str, delimiter=',', quotechar='"',
                 quoting=csv.QUOTE_MINIMAL):
    """Create a comma separated file from a table.

    After creating a table of model outputs using modelList_to_table
    this function can be used to create a CSV file from the table.

    Parameters
    table (list of lists): table to create CSV from
    csvfile (str): name of the CSV file to create
    delimiter (str): character to delimit CSV fields
    quotechar (str): character to use to quote field strings
    quoting (enum): quoting style to use (see csv.writer in Python docs)
    """
    with open(csvfile, 'w', newline='') as csvfile:
        out = csv.writer(csvfile, delimiter=delimiter,
                         quotechar=quotechar, quoting=quoting)
        for row in table:
            out.writerow(row)


def series_to_table(modelListSeries: ModelListSeries, header=True,
                    concat_names=None) -> List[List]:
    """Create a flat table from a series of model lists.

    This function "flattens" a time series of model lists so that each row in
    the table corresponds to a group in each model.  It's especially useful for
    analysing the model output, for example with numpy. It is also an interim
    step to converting a list of models to a comma-separated file.

    Example:
    table = series_to_table(simulate([my_model])

    Parameters
    modelListSeries (ModelListSeries): time series of model lists to convert
    header (bool): whether or not the first row should be a header
    concat_names (str): if not None then group names are concatenated,
                        separated by this string
    """
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
    """Create a comma separated file from a time series of model lists.

    This function is typically used in conjunction with "simulate"
    or "simulate_series" to generate a CSV file. E.g.

    series_to_csv(simulate(my_model), "mycsvfile.csv")

    Parameters
    modelListSeries (ModelListSeries): a time series of model lists
    csvfile (str): name of the CSV file to create
    header (bool): whether the first row of the csv file should be a header
    delimiter (str): character to delimit CSV fields
    quotechar (str): character to use to quote field strings
    quoting (enum): quoting style to use (see csv.writer in Python docs)
    concat_names (str): if not None then group names are concatenated,
                        separated by this string
    """
    table = series_to_table(modelListSeries, header, concat_names)
    table_to_csv(table, csvfile, delimiter, quotechar, quoting)


def simulate(modelList: ModelList, ident=None) -> ModelListSeries:
    """Iterate list of models and return a time series of model lists.

    Note that often the first parameter will only contain one model. It's
    typically only if you wish to iterate through distinct models that might be
    connected before or after each iteration through a hook in the
    "before_funcs" or "after_funcs" parameters that there would be more than
    one model in the list.

    Parameters:
    modelList (modelList): list of related models to iterate
    ident (int): unique identifier to use to identify this time series
                 (useful for generating a table or CSV file with multiple
                 model time series).

    """
    results = []
    for model in modelList:
        m = deepcopy(model)
        m['parameters'] = _make_parameters(m.get('parameters', {}))
        results.append(m)
    return _iterate_model(results, ident)


def _simulate(m):
    return simulate(m[0], m[1])


def simulate_series(modelListSeries: ModelListSeries,
                    processes=os.cpu_count()) -> ModelListSeries:
    """Execute series of models and return a time series of model lists.

    This function is useful for sensitivity analysis or calibration.
    The point of it is to execute the same or related models many times.
    It uses Python's multiprocessing library to execute the models in
    parallel. It returns a time series of model lists each with a unique
    identifier so that the output can be sorted appropriately afterwards.

    Parameters:
    modelListSeries (modelListSeries): series of model lists to execute in
                                       parallel
    processes (int): number of CPU processes to use (default uses one
                     process for each CPU on the machine)
    """
    mls_with_ident = [(m[0], m[1]) for m in
                      zip(modelListSeries, range(len(modelListSeries)))]
    with Pool(processes=processes) as pool:
        output = pool.map(_simulate, mls_with_ident)
    results = []
    for r in output:
        results += r
    return results
