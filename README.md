![georb_logo_v7](https://user-images.githubusercontent.com/102968112/182955140-320d8bae-7220-4e0c-8c9d-60481c1ce031.png)


# 

# Gravity and precisE ORBit determination system
  
GEORB is a software package for Precise Orbit Determination (POD) of Low Earth Orbiters (LEOs), gravity field recovery based on satellite gravity missions and design of future space missions. 

GEORB has been released as open source in 2022 while it has been developed by Thomas since 2007. 

GEORB has been used in orbit determination of the Gravity Recovery And Climate Experiment (GRACE), GRACE Follow-On (GRACE-FO) and Gravity Field and Steady-State Ocean Circulation (GOCE) missions. The current released version focuses on the POD and accelerometer calibration of the GRACE Follow-On mission.
 

---
## Cite GEORB: 

Thomas Loudis Papanikolaou (2023). GEORB: Release for precise orbit determination of low Earth orbiters and satellite gravity missions, *Software Impacts*, doi: [https://doi.org/10.1016/j.simpa.2023.100502](https://www.sciencedirect.com/science/article/pii/S2665963823000398). 


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


