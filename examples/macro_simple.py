from ziggy import macro

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
