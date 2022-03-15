# Catch-comparison-data
Data and code files for Detection and quantification of differences in catch rates among research vessel gears and commercial vessels

#### 

'Size curve analysis online.R' contains the statistical analyses of the catch comparison rates of the CPUE of scallops in 5 mm size groups. 

#### DATA FILES

'Comparison data for models.csv' is read in to the three analysis scripts.  

Length: shell width of king scallop in 5 mm groups. Groups are rounded down, so 110 group contains scallops 110 to 114 in shell width.  

Station: station name.  

Test: Number of king scallops caught by the designated ‘test’ gear-vessel. See ‘comp2’.  

Control: Number of king scallops caught by the designated ‘control’ gear-vessel. See ‘comp2’. 

Areatest: swept area of test in km2 

Areacont: swept area of control in km2 

Day: day of experiment where 1 to 4 correspondS to 25th to 28th April 2021.  

Depth: water depth in m as measured from the research vessel, which has the echosounder 3 m lower than sea level.  

Samp_test: raising factor for test. 

Samp_cont: raising factor for control.  

Cpuetest: test/areatest 

Cpuecont: control/areacont 

Cpuetotal: cpuetest + cpuecont 

Prop: cpuetest/cpuetotal 

Qratio: samp_test/samp_cont 

Comp2: name of comparison, with the name listed before the 'vs' the vessel designated as 'test' and the one after designated as 'control'

###################
###################

'Model data 3 knots.csv' contains size-structured results after model fitting.  

Sl: Length scaled and centred (see R ‘scale’ function for details).  

CpuetotalN, Station, day, qratio: dummy columns that can be ignored.  

Prop: predicted proportion of mature king scallops.  

Lo: 5% quantile of bootstrapped prop.  

Mid: 50% quantile of bootstrapped prop.  

Hi: 95% quantile of bootstrapped prop.  

Length: king scallop shell width in 5 mm classes. See above for details.  

Comp2: as above.  

Mod: Which model was selected by AIC? Currently, mod is a third order polynomial and mod2 is a basis spline with three knots.  

Ratio: predicted catch ratio of king scallops.  Not used in paper.

Cr_lo: 5% quantile of bootstrapped ratio.  Not used in paper.

Cr_mid: 50% quantile of bootstrapped ratio.  Not used in paper.

Cr_hi: 95% quantile of bootstrapped ratio. Not used in paper. 

##################
##################
'Model residuals 3 knots.csv' contains size-structured residual estimates from model fit. The majority of columns are defined in the previous datasets.  

sdepth: ignore.

Res: deviance residuals. 

############
############

'Model parameters 3 knots.csv' contains estimates of the model parameters for each comparison 

Intercept: model intercept
b1, b2, b3: beta 1:3 as described in paper

############
############

'Distance between vessels.csv" contains estimates of the distances between pairs of vessels during the hauls

dista: distance between vessels at starting coordinates in m 
station: haul name 
comp: comparsion

############
############
