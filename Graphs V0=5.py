# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 11:02:37 2023

@author: 61439
"""


import matplotlib.pyplot as plt 

#V0=5m/s
# thrust forces 

# line 1 points 
x1 = [2.8,11.0027383,16.871021,22.9641225,32.3076383,41.5730197,50.4110617,58.5344352,65.7500952,71.9674921,77.1859743,78.7133469,80.1402166,82.7084849,84.9251565,86.8264859,88.4486629,89.166] 

y1 = [39.08822781,87.69891405,155.3481109,574.6113278,963.8223107,1180.717498,1470.475758,1756.919243,1990.23276,2143.303317,2222.392353,2233.795108,2235.729066,2194.4397,2070.650846,1825.330288,1312.719673,0] 
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.', linestyle="dashed") 
  
# line 2 points 
x3 = [0,	2.64302,	5.37977,	8.20274,	11.1031,	14.071,	17.0953,	20.1641,	23.2647,	26.3837,	29.5076,	32.6228,	35.7156,	38.773,	41.7824,	44.732,	47.6111,	50.4099,	53.1201,	55.7344,	58.247,	60.6534,	62.9501,	65.1352,	67.2076,	69.1675,	71.0159,	72.7545,	74.386,	75.9133,	77.3402,	78.6705,	79.9085,	81.0585,	82.1252,	83.113,	84.0265,	84.8703,	85.6487,	86.366] 
x2=[value + 2.8 for value in x3]

y2 = [38.7712,	46.9196,	111.002,	124.794,	135.616,	155.27,	179.439,	383.207,	816.356,	925.107,	910.575,	999.703,	1087.22,	1179.13,	1272.64,	1367.93,	1466.71,	1565.16,	1659.34,	1751.97,	1838.59,	1916.88,	1985.36,	2043.1,	2094.57,	2135.39,	2172.95,	2196.96,	2218.59,	2227.71,	2228.9,	2216.76,	2189.06,	2138.68,	2063.8,	1967.6,	1824.93,	1618.68,	1306.95,	0] 
# plotting the line 2 points  
plt.plot(x2, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('Span [m]') 
# naming the y axis 
plt.ylabel('Thrust Force [N/m]') 
# giving a title to my graph 
plt.title('Thrust Force Along Blade - V0=5 m/s') 
  
# show a legend on the plot 
plt.legend() 
  
plt.savefig("Thrust-V0=5.pdf", format="pdf", bbox_inches="tight")
# function to show the plot 
plt.show() 

#Torque forces 
  
# line 1 points 
y1 = [-9.83145422605571,-26.5021710807464,	15.904542162128,	156.392505101153,	193.43864930329,	193.36835322703,	193.201427161578,	192.372039884353,	190.515863085074,	187.066036053448,	179.74417924332,	175.882795269235,	171.08864893269,	158.019105117614,	138.882084407882,	110.739961260839,	63.8221053814183,	0]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.', linestyle="dashed") 
  
# line 2 points 
y2 = [-9.90107,	-23.1806,	46.5055,	36.9497,	25.9886,	19.5773,	10.4662,	64.8713,	192.205,	191.533,	194.733,	194.171,	193.351,	195.708,	195.128,	194.622,	194.224,	193.788,	193.29,	192.739,	192.135,	191.475,	190.712,	189.782,	188.623,	187.154,	185.267,	182.844,	179.776,	175.919,	171.141,	165.275,	158.132,	149.487,	139.072,	126.465,	110.995,	91.377,	64.1308,	0]
# plotting the line 2 points  
plt.plot(x2, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('Span [m]') 
# naming the y axis 
plt.ylabel('Torque Force [N/m]') 
# giving a title to my graph 
plt.title('Torque Force Along Blade - V0=5 m/s') 
  
# show a legend on the plot 
plt.legend() 

plt.savefig("Torque-V0=5.pdf", format="pdf", bbox_inches="tight")
  
# function to show the plot 
plt.show() 