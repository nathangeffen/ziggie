# Ziggie macro modelling

THIS IS A WORK IN PROGRESS - NOT YET READY FOR USE

Ziggie is a Python package for infectious disease modelling.

The macro module facilitates compartmental modelling using difference equations,
or macro models.

## Quick start

Ziggie requires Python 3.5 or later.

First install it:

```bash
pip install ziggie
```

Then test that it's working:

```bash
python test.py
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

# Make a list of outputs with the results
# which by default are printed out every 50 days and
# at the beginning and end of the simulation.
# Each entry in the results table is an updated model
# For a particular day.
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

## A sophisticated example

Let's say we want to model the Covid-19 epidemic in South Africa.

Our model world looks like this:

- We have three distinct geographies: (1) urban areas with formal housing (2)
  urban areas with informal housing and (3) rural areas.
- Each area has three distinct age groups: (1) 0-24 years old, (2) 25-54, and
  (3) 55 and older.

TO DO
