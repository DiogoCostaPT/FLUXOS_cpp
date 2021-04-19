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

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).

Example:

.. code-block:: json

    {
        "PROJECT_NAME": "Test",
        "GEOGRAPHICAL_LOCATION": "Svalbard",
        "AUTHORS": "Diogo Costa",
        "DATE": "May_2020",
        "COMMENT": "Default file for testing purposes",
        "OPENWQ_INPUT": {
            "CONFIG_FILEPATH": "bin/openWQ_config.json",
            "BGC_CYCLES_FILEPATH": "bin/openWQ_BGC_cycling.json",
            "SINKSOURCE_FILES": {
                "1": {
                    "LABEL": "fertilizer_N",
                    "FILEPATH": "bin/openWQ_source_fertN.json"
                },
                "2": {
                    "LABEL": "fertilizer_P",
                    "FILEPATH": "bin/openWQ_source_fertP.json"
                }
            }
        },
        "openWQ_OUTPUT": {
            "RESULTS_FOLDERPATH": "bin/Output_OpenWQ",
            "FORMAT": "HDF5",
            "CHEMICAL_SPECIES":["SOIL_NO3","SOIL_NO3"],
            "COMPARTMENTS": ["SNOW","SOIL"],
            "TIMESTEP": [1,"h"]
        }
    }

