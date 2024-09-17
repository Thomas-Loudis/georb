![GEORB_logo_bbw_2024](https://github.com/Thomas-Loudis/georb/assets/102968112/17c74208-7d07-4549-9024-de82b546dc69)

# 

# Gravity and precisE ORBit determination system

[GEORB](https://georb.gr) is a software for Precise Orbit Determination (POD) of Low Earth Orbiters (LEOs) and satellite gravity missions, gravity field recovery and orbit design of future space missions. 

GEORB was created by Dr. [Thomas Loudis Papanikolaou](https://thomaspap.com) as a long-term project, initiated in 2007 and released as open source in 2022. It is designed as a versatile tool to support academic research in satellite geodesy and orbital mechanics, as well as industrial applications in Earth observation and space exploration. GEORB is currently exploited and further developed at [DORUS Space Lab](https://dorus.space). It has also been deployed at Aristotle University of Thessaloniki (GR) and Aalborg University Copenhagen (DK) and has supported projects of the European Space Agency (ESA) and Geoscience Australia.  

GEORB specialises in data analysis of satellite gravity missions and simulation of future space gravity missions, supporting revolutionary missions such as the [Gravity Recovery And Climate Experiment (GRACE)](https://www.jpl.nasa.gov/missions/gravity-recovery-and-climate-experiment-grace), and [Gravity Field and Steady-State Ocean Circulation (GOCE)](https://earth.esa.int/eogateway/missions/goce) missions. The current release delivers data products of precise orbits, calibrated accelerometer data and time-variable gravity field models for the NASA/GFZ’ [GRACE Follow-On](https://www.jpl.nasa.gov/missions/gravity-recovery-and-climate-experiment-follow-on-grace-fo) mission. 

Data Products (major):
- Precise Orbit data
- Gravity Field models
- Calibration of on-board accelerometer data
- Orbit simulation data
 

---
## Cite GEORB: 

Thomas Loudis Papanikolaou (2023). GEORB: Release for precise orbit determination of low Earth orbiters and satellite gravity missions, *Software Impacts*, doi: [https://doi.org/10.1016/j.simpa.2023.100502](https://www.sciencedirect.com/science/article/pii/S2665963823000398). 

---
## Data Products: 

Thomas Loudis Papanikolaou (2023b). Precise orbit data and Gravity field models for GRACE Follow-On mission. VBN. doi: [https://doi.org/10.5278/8b4d5e8d-0010-4461-a58e-b5b39637bd51](https://doi.org/10.5278/8b4d5e8d-0010-4461-a58e-b5b39637bd51)

Thomas Loudis Papanikolaou (2023c). DORUS time-variable gravity field models. VBN. doi: [https://doi.org/10.5278/b9644d2d-697f-4b73-bf73-c79582e6572b](https://doi.org/10.5278/b9644d2d-697f-4b73-bf73-c79582e6572b)


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

2. Download the models' data required by executing the script file `georb_data_models.m` stored in the folder `'../scripts/'`. The data will be save in the folder `'../data/'`.

```
cd scripts/
georb_data_models
```

3. Satellite Missions data are required when operating in the `'orbit missions mode'` and other modes requiring satellite data.


The current release supports the GRACE missions.

Data of the GRACE Follow-On and GRACE missions can be accessed from the official data centers servers as follows: 

JPL/NASA: https://podaac.jpl.nasa.gov/GRACE-FO

GFZ: https://isdc.gfz-potsdam.de/grace-fo-isdc/

The GRACE/GRACE-FO data need to be stored in the folder `'../data/'` 
 

4. Execute the GEORB main script `georb_main.m` in the folder `'../main/'` 

```
cd main/
georb_main
```

5. The computed restuls are written to output data (ascii) files saved in the folder `'../data_output/'`


# 


---


