# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 11:39:58 2023

@author: 61439
"""


import matplotlib.pyplot as plt 

#V0=9m/s
# thrust forces 

# line 1 points 
x1 = [2.8,11.0027383,16.871021,22.9641225,32.3076383,41.5730197,50.4110617,58.5344352,65.7500952,71.9674921,77.1859743,78.7133469,80.1402166,82.7084849,84.9251565,86.8264859,88.4486629,89.166] 

y1 = [126.645858106526,	284.144481507937,	503.327879391244,	1861.74070203203,	3122.78428662663,	3825.5246943402,	4764.34145726759,	5692.41834670399,	6448.35414101454,	6944.30274732172,	7200.55122451767,	7237.49615013635,	7243.76217262896,	7109.98462685754,	6708.90873969693,	5914.07013411167,	4253.21174068392,	0] 
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.', linestyle="dashed") 
  
# line 2 points 
x3 = [0,	2.64302,	5.37977,	8.20274,	11.1031,	14.071,	17.0953,	20.1641,	23.2647,	26.3837,	29.5076,	32.6228,	35.7156,	38.773,	41.7824,	44.732,	47.6111,	50.4099,	53.1201,	55.7344,	58.247,	60.6534,	62.9501,	65.1352,	67.2076,	69.1675,	71.0159,	72.7545,	74.386,	75.9133,	77.3402,	78.6705,	79.9085,	81.0585,	82.1252,	83.113,	84.0265,	84.8703,	85.6487,	86.366] 
x2=[value + 2.8 for value in x3]

y2 = [125.668,	152.203,	359.984,	404.142,	439.09,	502.852,	580.416,	1245.87,	2645.95,	2998.12,	2950.81,	3239.5,	3523.09,	3821.15,	4124.36,	4433.4,	4753.7,	5072.95,	5378.35,	5678.73,	5959.64,	6213.58,	6435.81,	6623.14,	6790.09,	6922.54,	7044.39,	7122.34,	7192.57,	7222.24,	7226.19,	7186.94,	7097.2,	6933.96,	6691.26,	6379.44,	5916.95,	5248.31,	4237.7,	0]
# plotting the line 2 points  
plt.plot(x2, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('Span [m]') 
# naming the y axis 
plt.ylabel('Thrust Force [N/m]') 
# giving a title to my graph 
plt.title('Thrust Force Along Blade - V0=9 m/s') 
  
# show a legend on the plot 
plt.legend() 

plt.savefig("Thrust-V0=9.pdf", format="pdf", bbox_inches="tight")

# function to show the plot 
plt.show() 

#Torque forces 
  
# line 1 points 
y1 = [-31.8539116924203,	-85.8670343016174,	51.5307166052952,	506.711716527736,	626.741223742662,	626.513464455582,	625.972624003516,	623.285409225308,	617.271396395644,	606.093956813174,	582.371140748361,	569.860256672325,	554.327222541918,	511.981900581071,	449.97795348154,	358.797474485119,	206.783621435796,	0]
# plotting the line 1 points  
plt.plot(x1, y1, label = "Python", marker ='.', linestyle="dashed") 
  
# line 2 points 
y2 = [-32.6927,	-75.9233,	149.669,	118.603,	83.4283,	62.8967,	32.8771,	211.372,	621.148,	619.038,	629.56,	627.795,	625.198,	632.919,	631.1,	629.523,	628.27,	626.884,	625.296,	623.531,	621.596,	619.484,	617.046,	614.069,	610.344,	605.621,	599.542,	591.733,	581.832,	569.376,	553.935,	534.971,	511.867,	483.902,	450.203,	409.4,	359.325,	295.819,	207.609,	0]
# plotting the line 2 points  
plt.plot(x2, y2, label = "Ashes", marker ='.', linestyle="dashed") 
  
# naming the x axis 
plt.xlabel('Span [m]') 
# naming the y axis 
plt.ylabel('Torque Force [N/m]') 
# giving a title to my graph 
plt.title('Torque Force Along Blade - V0=9 m/s') 
  
# show a legend on the plot 
plt.legend() 

plt.savefig("Torque-V0=9.pdf", format="pdf", bbox_inches="tight")

  
# function to show the plot 
plt.show() 