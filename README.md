# Ziggie macro modelling

Ziggie is a Python package for infectious disease modelling.

The macro module facilitates compartmental modelling using difference equations
(i.e. macro models).

## Quick start

Ziggie requires Python 3.5 or later.

First install it:

```bash
pip install ziggie
```

Here's Python code to create and run a simple SIR model

```Python
from ziggie import macro

# Simple SIR model with susceptible population of 1 million and
# one infection. Effective contact rate per day is 0.6 and the
# infection duration is ten days.

simple = {
    'name': 'Simple model',
    'compartments': {
        'S': 1000000, # Susceptible
        'I': 1,       # Infectious
        'R': 0        # Recovered
    },
    'transitions': {
        # (sometimes called beta in the literature)
        'S_I': 0.6, # Effective contact rate
        'I_R': 0.1  # Recovery rate per day
    },
}

""" Make a list of outputs with the results
    which by default are stored every 50 days and
    at the beginning and end of the simulation.
    Each entry in the results table is an updated model
    for a particular iteration (which corresponds to a day).
"""
results = macro.simulate([simple])
print("Final day's results")
print(results[-1])

# Flatten the final day's results into a table
table = macro.modelList_to_table(results[-1])
print("Time series table")
print(table)

# Put all the results into a CSV file
macro.series_to_csv(results, "my_csv_file.csv")

# Run model and create CSV in one step
macro.series_to_csv(macro.simulate([simple]), "my_csv_file.csv")

```

The output is:

```Python
Final day's results
[{'name': 'Simple model', 'compartments': {'S': 65525.67886409458, 'I': 1.5383929100326267e-07, 'R': 56934475.32113576}, 'transitions': {'S_I': 0.6, 'I_R': 0.1}, 'parameters': {'from': 0, 'to': 365, 'record_frequency': 50, 'reduce_infectivity': 1.0, 'asymptomatic_infectiousness': 1.0, 'treatment_infectiousness': 1.0, 'noise': 0.0, 'discrete': False, 'record_first': True, 'record_last': True, 'transition_funcs': {'S_I': <function delta_S_I at 0x7fa842cdb8b0>, 'S_E': <function delta_S_I at 0x7fa842cdb8b0>, 'S_I1': <function delta_S_I1 at 0x7fa84651e5e0>, 'S_E1': <function delta_S_I1 at 0x7fa84651e5e0>, 'B_S': <function delta_birth_X at 0x7fa84651e3a0>, 'default': <function delta_X_Y at 0x7fa846519310>}, 'before_funcs': [], 'after_funcs': []}, 'iteration': 365}]
Time series table
[['iter', 'name_0', 'S', 'I', 'R'], [0, 'Simple model', 57000000, 1, 0], [50, 'Simple model', 2200495.449318898, 28345727.9672264, 26453777.583454713], [100, 'Simple model', 66701.56917131442, 167716.11455651774, 56765583.316272154], [150, 'Simple model', 65531.91780545574, 898.2457754223102, 56933570.83641911], [200, 'Simple model', 65525.71227208052, 4.810125157065176, 56934470.477602795], [250, 'Simple model', 65525.67904299378, 0.025758303410947692, 56934475.29519872], [300, 'Simple model', 65525.67886505151, 0.00013793615973905623, 56934475.32099703], [350, 'Simple model', 65525.67886409864, 7.38650518227673e-07, 56934475.32113517], [365, 'Simple model', 65525.67886409458, 1.5383929100326267e-07, 56934475.32113576]]
```

You can also open one of the generated CSV files in your favourite spreadsheet program.

It's also easy to use with numpy. In the above example the table variable
consists of string and floats. Numpy arrays really only make sense if they're
numbers, so we'll chop off the strings before conversion.. Here's a continuation
of the above example:

```Python
import numpy

table_float = [t[2:] for t in table[1:]]
print(numpy.array(table_float))

```

The output should be:

```Python
array([[5.70000000e+07, 1.00000000e+00, 0.00000000e+00],
       [2.20049545e+06, 2.83457280e+07, 2.64537776e+07],
       [6.67015692e+04, 1.67716115e+05, 5.67655833e+07],
       [6.55319178e+04, 8.98245775e+02, 5.69335708e+07],
       [6.55257123e+04, 4.81012516e+00, 5.69344705e+07],
       [6.55256790e+04, 2.57583034e-02, 5.69344753e+07],
       [6.55256789e+04, 1.37936160e-04, 5.69344753e+07],
       [6.55256789e+04, 7.38650518e-07, 5.69344753e+07],
       [6.55256789e+04, 1.53839291e-07, 5.69344753e+07]])
```

## Transitions

A transition name consists of two compartment names separated by an
underscore. There are default transition functions, depending on the
compartments involved, which you can override by including a 'parameters'
dictionary in your model, and then overriding entries in its 'transition_funcs' sub-dictionary.

This is what the default *transition_funcs* dictionary looks like:

```Python
    'transition_funcs': {
        'S_I': delta_S_I,
        'S_E': delta_S_I,
        'S_I1': delta_S_I1,
        'S_E1': delta_S_I1,
        'B_S': delta_birth_X,
        'default': delta_X_Y
    }
```

The transition function to execute is looked up in this dictionary. If it's not
found the function specified by 'default' is executed.

For most transitions, a proportion of the "from" compartment is moved
to the "to" compartment. The delta_X_Y function takes care of this.

But for new infections you almost always want a standard SIR-like equation such
as this:

delta = susceptibles * number of contacts per iteration *
                       risk of infection per contact *
                       total number of infectious individuals /
                       total population

susceptibles -= delta
infectious (or exposed) individuals += delta

Two functions are provided to deal with this: the very simple delta_S_I
and the more sophisticated (but slower) delta_S_I1.

The delta_SI1 function calculates the total number of infectious individuals by
adding all the compartments starting with an *I* as well as all compartments
starting with an *A* (asymptomatic individuals) and *T* (treated
individuals). Moreover the infectiousness of the asymptomatic individuals is
multiplied by the parameter *asymptomatic_infectiousness* and the infectiousness
of treated individuals by *treatment_infectiousness* (both default to 1).

## Compartment names which have meaning

Some of the compartment name prefixes are meaningful, in that the code might
make assumptions about the compartment. A compartment name generally starts with
one of these meaningful prefixes and then a unique identifier. E.g. I1, I2, I3
or I4 for various stages of infectiousness.

- S - Susceptible  (See delta_S_I and delta_S_I1)
- E - Exposed  (See delta_S_I and delta_S_I1)
- I - Infectious (See delta_S_I and delta_S_I1)
- A - Asymptomatic (See delta_S_I1 and asymptomatic_infectiousness)
- T - On treatment (See delta_S_I1 and treatment_infectiousness)
- N - Population size. Strictly reserved. Do not prefix any compartment N.
- D - Dead (not included in totalling N)
- B - Birth (not included in totalling N)
- R - Recovered
- M - Maternal immunity
- V - Vaccinated

You can also prefix a compartment name with any other letter.

## Noise

You can also add noise to your model so that it is stochastic. We add a
'parameters' key to our model, and within the parameters, we add a 'noise'
key. Like so:

```Python
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
    'parameters': {
        'noise': 0.05
    }
}
```

Now every transition calculation is multipled by a uniform random number in the
range [1-0.05, 1+0.05].

## Parameters

Besides 'noise' there are many other parameters that can be modified
directly. These are the default parameters:

```Python
{
    # The model is executed from iteration 'from' to iteration 'to' - 1.
    'from': 0,
    'to': 365,
    # The results are recorded every 'record_frequency' iterations. Set
    # to 1 if you want to record the output of every iteration.
    'record_frequency': 50,  # Write results every X iterations
    # Multiply S_E or S_I by X every iteration if reduce_infectivity
    # function executed. This is useful for modelling heterogeneity, i.e. the
    # fact that usually the most susceptible people get infected earliest in an
    # epidemic
    'reduce_infectivity': 1.0,
    # If you have multiple infectiousness compartments including
    # asymptomatic and treatment compartments, you can
    'asymptomatic_infectiousness': 1.0,
    'treatment_infectiousness': 1.0,
    # Add stochastic noise to transitions
    'noise': 0.0,
    # Round all transition calculations to round numbers (not tested yet)
    'discrete': False,
    # Include the initial state of the model in the results
    'record_first': True,
    # Include the final state of the model in the results.
    # If record_frequency is 1 or divides into the "to" parameter,
    # you probably want to set this to False.
    'record_last': True,
    # The transition functions
    'transition_funcs': {
        'S_I': delta_S_I,
        'S_E': delta_S_I,
        'S_I1': delta_S_I1,
        'S_E1': delta_S_I1,
        'B_S': delta_birth_X,
        'default': delta_X_Y
    },
    # Any functions specified here are executed for each model before
    # each iteration
    'before_funcs': [],
    # Any functions specified here are executed for each model after
    # each iteration
    'after_funcs': [],
}
```

## More sophisticated example

Let's say we want to model the Covid-19 epidemic in South Africa.

Our model world might like this:

- We have three distinct geographies: (1) urban areas with formal housing (2)
  urban areas with informal housing and (3) rural areas.
- Each area has three distinct age groups: (1) 0-24 years old, (2) 25-54, and
  (3) 55 and older.
- People start off (S)usceptible. Upon being infected they are (E)xposed. They
  transition from (E)xposed to either (A)symptomatic or infectious with mild
  symptoms(Im). Asymptomatic people transition to (R)ecovered. Infectious with
  mild symptoms transition to infectious with critical symptoms (Ic) and then
  either to (R)ecovered or (D)eath.
- After each iteration, there is a bit of migration between the three
  geographies.

The source code file ziggie/samples.py contains a class called MacroModels
which contains a method called corona that implements this. Here it is:

```Python
from copy import deepcopy
from ziggie import macro

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

m = corona()
results = macro.simulate(m)
print("Number of results:", len(results)) # Outputs 366
print("Number of models:", len(results[-1])) # Outputs 3
totals = [macro.calc_totals(results[-1][i])
          for i in range(len(results[-1]))
print(macro.sum_totals(totals))
```

The output is something like this:

```
Number of results: 366
Number of models: 3
{'N': 59685873.7194506, 'S': 47387513.37539025, 'E': 283059.79850507394, 'Im': 108802.85166144818, 'Ic': 106605.46023480814, 'A': 264806.9195997677, 'R': 11535085.314059254, 'D': 134156.28054939664}
```

# Reference

## Main data structures

A model specification consists of groups. A model is itself the
top-level group.

Each group can itself contain groups.

A group specification is a dictionary containing or more of the following
fields:

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


A ModelList is a list of models. The *simulate* function takes
a list of models as its first parameter. Why? Let's say you have
distinct regions you'd like to model, say London and Manchester.
You would specify a model for each and then put them in a list (i.e. you
now have a ModelList data structure). You can then run "simulate"
on this ModelList.

The output is two distinct sets of independent new model specification
with updated compartments on each iteration.  i.e. a list of ModelLists.

This is called a ModelListSeries.

To summarise from lowest level to highest level:

- Group: A dictionary defining part of a model. A model is itself a group.

- Model: The highest level group. It defines a complete model.

- ModelList: A list of models. This is what you pass to the *simulate* function

- ModelListSeries: A list of ModelLists. This is the output of *simulate*.


## Function reference

### Most frequently used functions

simulate(modelList: ModelList, ident=None) -> ModelListSeries

    Iterate list of models and return a time series of model lists.

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

simulate_series(modelListSeries: ModelListSeries, processes=int)
                -> ModelListSeries

    Execute series of models and return a time series of model lists.

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


series_to_csv(modelListSeries: ModelListSeries,
              csvfile: str, header=True, delimiter=',',
              quotechar='"', quoting=0, concat_names=None)

    Create a comma separated file from a time series of model lists.

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

series_to_table(modelListSeries: ModelListSeries, header=True,
                concat_names=None) -> List[List]

    Create a flat table from a series of model lists.

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



calc_totals(model: Model) -> Dict[str, float]

    Calculate sum of each compartment in all groups and return dict.

    This function traverses a model and calculates the sum of each compartment
    across all the groups. It also calculates N, the total living population.

grand_sum_totals(totals: List[Dict[str, float]], ignore=['B', 'N']) -> float

    Calculate the sum of all compartments.

    Takes the output list generated using calc_totals and calculates the sum of
    all the compartments.

    E.g. print(grand_sum_totals([calc_totals[m] for m in model_list])

    Parameters:

    totals (output of calc_totals): list of compartment totals across models
    ignore (list of str): list of compartments to ignore


### Supporting functions

modelList_to_table(modelList: List[Dict], header=True,
                   concat_names=None) -> List[List]

    Create and return a table from a list of models.

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

model_to_table(model: Dict, concat_names=None) -> List[List]

    Create and return a table from a model.

    This converts a model to list of lists, where each list corresponds
    to a group in the model. This is also an interim step to converting
    a model output to a csv file.

    Parameters
    model (Model): model to convert
    concat_names (str): if not None then group names are concatenated,
                        separated by this string

reduce_infectivity(model: Dict, modelList=None)

    Reduce the values of the effective contact rates of a model.

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


sum_infectiousness(model: Dict) -> float

    Return weighted sum of all infection compartments.

    The infection compartments are all those that begin with an 'I'
    (infectious), 'A' (asymptomatic) or 'T' (treatment).
    Each 'A' compartment is multipled by the parameter
    asymptomatic_infectiousness. Each 'T' compartment is multiplied by
    the parameter treatment_infectiousness.

    Parameters:

    model: the model to sum

sum_totals(totals: List[Dict[str, float]]) -> List[Dict[str, float]]

    Calculate the sum of compartments across multiple models.

    Takes the output list generated using calc_totals and calculates the sum of
    each the compartments.

    E.g. print(sum_totals([calc_totals[m] for m in model_list])

    Parameters:

    totals (output of calc_totals): list of compartment totals across models

table_to_csv(table: List[List], csvfile: str,
             delimiter=',', quotechar='"', quoting=0)

    Create a comma separated file from a table.

    After creating a table of model outputs using modelList_to_table
    this function can be used to create a CSV file from the table.

    Parameters

    table (list of lists): table to create CSV from
    csvfile (str): name of the CSV file to create
    delimiter (str): character to delimit CSV fields
    quotechar (str): character to use to quote field strings
    quoting (enum): quoting style to use (see csv.writer in Python docs)

traverse(group: Dict) -> Generator[Dict, NoneType, NoneType]

    Generate all the groups of a model or group.

    Parameters:

    group: the group whose subgroups to traverse

### Transition functions

These are the provided transition functions, but you can easily code your
own and modify the *transition_funcs* dictionary in *parameters* to make
the system aware of it.

All have the same function signature:

    transition - usually named from_to
    float  - this is taken from "transitions" entry of a model group.
    compartments - compartments dictionary for a model group
    totals - dictionary of the model level sum of all the compartments in this
             group, including N the population size)
    model - the entire current state of the model


delta_S_I(from_to, beta, compartments, totals, model=None)

    Return number of new infections.

    Parameters:

    from_to (str): transition name consisting of two compartment names
                    separated by an underscore (e.g. S_Ic)
    beta (float): effective contact rate
    compartments (dict): dictionary of compartments including the two
                         specified in from_to
    totals (dict): dictionary containing a key 'N' that is the sum of the
                   total population for this model.
    model (Model): Unused but part of function signature

delta_S_I1(from_to, beta, compartments, totals, model=None)

    Return number of new infections.

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

delta_X_Y(from_to, prop, compartments, totals, model=None)

    Return number individuals to be moved from one compartment to another.

    Parameters:

    from_to (str): transition name consisting of two compartment names
                    separated by an underscore (e.g. I_R)
    prop (float): proportion of "from" compartment to move
    compartments (dict): dictionary of compartments including the two
                         specified in from_to
    totals (dict): Unused but part of function signature
    model (Model): Unused but part of function signature

delta_birth_X(from_to, prop, compartments, totals, model=None)

    Return number new individuals in a population.

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

### Experimental

R0(modelListSeries: List[List[Dict]]) -> List[float]

    Calculate R0 for a time series of model outputs.

    Highly inaccurate and needs to be reconceptualised. Don't use for now.
