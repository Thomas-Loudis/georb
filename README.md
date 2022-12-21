![georb_logo_v7](https://user-images.githubusercontent.com/102968112/182955140-320d8bae-7220-4e0c-8c9d-60481c1ce031.png)


# 

# Gravity and precisE ORBit determination system
  
GEORB is a software for Precise Orbit Determination (POD) of Low Earth Orbiters (LEOs), data analysis of satellite gravity missions and design of future space missions. 	
It has been mostly applied in precise orbit determination of the Gravity Recovery And Climate Experiment (GRACE), GRACE Follow-On (GRACE-FO) and Gravity Field and Steady-State Ocean Circulation (GOCE) missions.   


---
## Cite GEORB: 

Papanikolaou, T. (2022). Precise orbit determination and accelerometer data modelling of the GRACE Follow-On mission, GRACE/GRACE-FO Science Team Meeting 2022, Potsdam, Germany, 18–20 Oct 2022, GSTM2022-90, https://doi.org/10.5194/gstm2022-90, 2022.

---
## Guide: Instructions for configuration and data requirements

GEORB can be executed through applying the following steps:

1. Set the configuration files located in the folder `'../config/'`. Description of the configurable parameters is provided within the config files.  

`main_config.in` :: Master Configuration file for setting the basic modes 

`orbit_model.in` :: Configuration file for the orbit modelling and methods

`ic_config.in`   :: Initial Consditions file for setting the list of the Satellites/Objects along with the Initial Epoch per object 

```
cd config/
edit main_config.in
cd ..
```

```
cd config/
edit orbit_model.in
cd ..
```

```
cd config/
edit ic_config.in
cd ..
```

2. Download the models' data required by executing the script file `georb_data_models.m` stored in the folder `'../scripts/'`

```
cd scripts/
georb_data_models
```

3. Satellite Missions data are required when operating in the `'orbit missions mode'` and other modes related to missions. 

The current release supports the GRACE missions.

Download the data of the GRACE Follow-On and GRACE missions from the official data centers servers: 

JPL/NASA: https://podaac.jpl.nasa.gov/GRACE-FO

GFZ: https://isdc.gfz-potsdam.de/grace-fo-isdc/


4. Execute the main script `georb_main.m` in the folder `'../main/'` 

```
cd main/
georb_main
```

5. The computed restuls are written to output data (ascii) files saved in the folder `'../data_output/'`


# 


---
## References:

Papanikolaou T. (2012). Dynamic modelling of satellite orbits in the frame of contemporary satellite geodesy missions, Ph.D. Dissertation, Aristotle University of Thessaloniki (AUTH), Greece.

Papanikolaou T., Tsoulis D. (2016). Assessment of numerical integration methods in the context of low Earth orbits and inter-satellite observation analysis, Acta Geodetica et Geophysica, Vol. 51, No. 4, pp. 619-641, https://doi.org/10.1007/s40328-016-0159-3.

Papanikolaou T., Tsoulis D. (2018). Assessment of Earth gravity field models in the medium to high frequency spectrum based on GRACE and GOCE dynamic orbit analysis, Geosciences, 8(12): 441 , In: Special Issue "Gravity Field Determination and Its Temporal Variation" (Eds. M. Bagherbandi), https://doi.org/10.3390/geosciences8120441.

Papanikolaou T. (2022). Precise Orbit Determination and accelerometry calibration modelling of the GRACE Follow-On mission, Nordic Geodetic Commission (NKG) General Assembly, 5-8 September 2022, Copenhagen, Denmark.
