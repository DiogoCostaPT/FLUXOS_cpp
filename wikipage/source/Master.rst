Master configuration file
==================================

The master configuration is JSON file that provides FLUXOS-OVERLAND with information and instructions to run a simulation. It contains pairs of variable/command names and values. The variables/commands needed to run a simulation are:

.. list-table:: Variables/Keywords in MASTER file (JSON type)

    *   - Variable/Keyword
        - any useful descriptive information about the project
    *   - DEM_FILE
        - full path to the ASCII-ERSI DEM file
    *   - PRINT_STEP
        - temporal resolution of the output files (this is not the model timestep as that is defined by the Courant–Friedrichs–Lewy Condition of numerical stability)
    *   - ROUGNESS_HEIGHT
        - average roughness of the terrain measured in height (m), it is used in friction model and water will be retained until this height
    *   - SOIL_RELEASE_RATE
        - release of soil contaminants (mg/hour)
    *   - SOIL_CONC_BACKGROUND
        - initial surficial soil concentration (mg/l)
    *   - SWE_STD
        - standard deviation of pre-melt SWE (for calculation of runoff-soil interaction during snowmelt, used the WINTRA model)
    *   - SWE_MAX
        - maximum SWE (for calculation of runoff-soil interaction during snowmelt, used the WINTRA model)
    *   - METEO_FILE (optional)
        - full path to the snowmelt/rainfall timeseries
    *   - INFLOW_FILE (optional)
        - full pa/rainfall timeseries

Example:

.. code-block:: json

    {
        "COMMNET": "Batch_1 - additional tests for paper",
        "DEM_FILE": "bin/Rosa_2m.asc",
        "INFLOW_FILE": {
            "FILENAME": "bin/Flow_forcing.fluxos",
            "DISCHARGE_LOCATION":{
                "NROW": 201,
                "NCOL": 310
            }
        },
        "PRINT_STEP": 3600,
        "ROUGNESS_HEIGHT": 0.005,
        "SOIL_RELEASE_RATE": 0.08,
        "SOIL_CONC_BACKGROUND": 1,
        "SWE_STD": 9.5,
        "SWE_MAX": 9.5
    }