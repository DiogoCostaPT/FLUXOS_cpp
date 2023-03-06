Master configuration file
==================================

The master configuration is JSON file that provides FLUXOS-OVERLAND with information and instructions to run a simulation. It contains pairs of variable/command names and values. The variables/commands needed to run a simulation are:

.. list-table:: Variables/Keywords in MASTER file (JSON type)

    *   - Variable/Keyword
        - any useful descriptive information about the project
    *   - ``DEM_FILE``
        - full path to the ASCII-ERSI DEM file
    *   - ``SIM_DATETIME_START``
        - Simulation start datetime
    *   - ``RESTART``
        - Restart option. If ``TRUE``, the model will load the last output file and restart from there.
    *   - ``PRINT_STEP``
        - temporal resolution of the output files (this is not the model timestep as that is defined by the Courant–Friedrichs–Lewy Condition of numerical stability)
    *   - ``ROUGNESS_HEIGHT``
        - average roughness of the terrain measured in height (m), it is used in friction model and water will be retained until this height
    *   - ``SOIL_RELEASE_RATE``
        - release of soil contaminants (mg/hour)
    *   - ``SOIL_CONC_BACKGROUND``
        - initial surficial soil concentration (mg/l)
    *   - ``SWE_STD``
        - standard deviation of pre-melt SWE (for calculation of runoff-soil interaction during snowmelt, used the WINTRA model)
    *   - ``SWE_MAX``
        - maximum SWE (for calculation of runoff-soil interaction during snowmelt, used the WINTRA model)
    *   - ``METEO_FILE`` (optional)
        - full path to the snowmelt/rainfall timeseries
    *   - ``INFLOW_FILE`` (optional)
        - full pa/rainfall timeseries
    *   - ``EXTERNAL_MODULES`` => ``ADE_TRANSPORT`` => ``STATUS``
        - activate 2D advection-dispersion solver
    *   - ``EXTERNAL_MODULES`` => ``ADE_TRANSPORT`` => ``D_COEF``
        - dispersion coefficient
    *   - ``EXTERNAL_MODULES`` => ``OPENWQ`` => ``STATUS``
        - activate openwq (disabled if ``ADE_TRANSPORT:STATUS`` is ``FALSE``)
    *   - ``EXTERNAL_MODULES`` => ``OPENWQ`` => ``MASTERFILE_DIR``
        - path to openwq master file
    *   - ``EXTERNAL_MODULES`` => ``WINTRA`` => ``STATUS``
        - activate wintra (disabled if ``ADE_TRANSPORT:STATUS`` is ``FALSE``)
    *   - ``EXTERNAL_MODULES`` => ``WINTRA`` => ``SWE_MAX``
        - max snow water equivalent for wintra calculations
    *   - ``EXTERNAL_MODULES`` => ``WINTRA`` => ``SWE_STD``
        - standard-deviation of snow water equivalent for wintra calculations

The JSON file supports C/C++ syntax for comments: single-line comment (``//``) or comment blocks (``/*`` and ``*/``).

Example:

.. code-block:: json

    {
        "COMMNET": "Batch_1 - additional tests for paper",
        "DEM_FILE": "Rosa_2m.asc",
        "SIM_DATETIME_START": "2017SEP01 12:15:00",
        "RESTART": false,

        "INFLOW_FILE": {
            "FILENAME": "Flow_forcing.csv",
            "DISCHARGE_LOCATION":{
                "X_COORDINATE": 425894,
                "Y_COORDINATE": 5785553
            }
        },

        "EXTERNAL_MODULES":{
            "ADE_TRANSPORT":{
                "STATUS": true,
                "D_COEF": 0.01
            },
            "WINTRA":{
                "STATUS": false,
                "SWE_STD": 9.5,
                "SWE_MAX": 9.5
            },
            "OPENWQ": {
                "STATUS": true,
                "MASTERFILE_DIR": "../openwq_in/openWQ_master.json"
            }
        },
        "METEO_FILE": "Qmelt_synthetic.fluxos",
        "OUTPUT": {
            "OUTPUT_FOLDER": "../fluxos_out",
            "PRINT_STEP": 3600,
            "H_MIN_TO_PRINT": 0.005
        },
        "ROUGNESS_HEIGHT": 0.005,
        "SOIL_RELEASE_RATE": 0.0,
        "SOIL_CONC_BACKGROUND": 0.0

    }



