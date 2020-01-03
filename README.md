# OceanExtremeEvents

Functions to be used for analysis of marine heatwaves and eddy-specific marine heatwaves.

# Contents

| File | Description |
| ---- | ----------- |
| eddySuppli_functions_G.py | Functions for eddy specific heatwave analysis |
| heatwave_functions_G.py | Functions for general heatwave analysis |
| hwPlot_functions_G.py | Heatwave data plotting|
| companion_G.py | File for setting parameter values, loading data and executing functions from above files |

# Additional Info 

**N.B.** functions in eddySuppli_functions.py build on tracked eddy data from Eric Oliver (https://github.com/ecjoliver/eddyTracking)

**N.B.** functions in heatwave_functions_G.py use some functions that are based on those by Eric Oliver (https://github.com/ecjoliver/marineHeatWaves)
   - Expands on Eric Oliver code by applying the methods to all grid cells within a given climatological region from ocean climate model output.
