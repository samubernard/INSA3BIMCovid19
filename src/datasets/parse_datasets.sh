#!/bin/bash

# ---------------
# Confirmed cases
# extract headers
head -1 time_series-ncov-Confirmed.csv  >confirmed_cases_metr_france.csv

# extract lines from Metropolitan France (France,France)
sed -n '/^\s*France\s*,\s*France/p' time_series-ncov-Confirmed.csv >>confirmed_cases_metr_france.csv

# ---------------
# Death cases  
# extract headers
head -1 time_series-ncov-Deaths.csv  >death_cases_metr_france.csv

# extract lines from Metropolitan France (France,France)
sed -n '/^\s*France\s*,\s*France/p' time_series-ncov-Deaths.csv >>death_cases_metr_france.csv

# ---------------
# Death cases  
# extract headers
head -1 time_series-ncov-Recovered.csv  >recovered_cases_metr_france.csv

# extract lines from Metropolitan France (France,France)
sed -n '/^\s*France\s*,\s*France/p' time_series-ncov-Recovered.csv >>recovered_cases_metr_france.csv

