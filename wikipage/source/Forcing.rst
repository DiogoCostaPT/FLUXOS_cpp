Forcing input file
==================================

    Metereological forcing to the model is provided through timeseries of rainfall and/or snowmelt. The units are mm/day, and the forcing will be applied uniformely across the entire model domain. The first row is for the header and is skipped.
    
    Example:
    
        t (sec), inflow/rain/snowmelt (mm/day), concentration (mg/l)

        0,    0.151, 2

        900,  0.144, 2

        1800, 0.158, 2

        2700, 0.181, 2

	... 
