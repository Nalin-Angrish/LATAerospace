''''

We try to fix the error that had arisen from the previous try by correcting the equations such that they work for a cambered airfoil
I will be focusing on the slope term first to see if the issue was caused by some error in that as I think that this is the only term
that is different in the cambered airfoils and the non cambered airfoils.

Here we are using some approximations such as-
1) Considering that the angles of attack are very small <<1 rad(around 57.2 degrees). Generally we do not increase the aoa more than 10 degrees
2)The airfoil is thin
3)We do not use the camer line for the vortex sheet, but instead, we use the x axis. Similar is the case for the normalcy condition
4)Small angle approximations are used for the free stream velocity
5)The airfoil is considered to be flat
'''

import csv
import matplotlib.pyplot as plt
import math
import numpy as np

def returnCamberLine(camberLineFileLocation):
    with open(camberLineFileLocation, mode='r', encoding='utf-8') as file:
        # Create a csv.reader object to iterate over lines in the CSV file
        csv_reader = csv.reader(file)
        header = next(csv_reader)
        xCamber = []
        yCamber = []

        x_index = header.index('X')
        y_index = header.index('Y')

        for row in csv_reader:
            xCamber.append(float(row[x_index]))
            yCamber.append(float(row[y_index]))

        
    
    
    return(xCamber, yCamber)

#Angle of Attack in radians
aoa = 10 * math.pi / 180

#Conditions
fsv = 20 #free stream velocity in m/s
rho = 1.225
cl = 0.0 #The coefficient of lift

#flap conditions
percFlap = 50 #This is the percentage of the wing that is the flap (measured from the trailing edge)
flapAngle = 5 * math.pi / 180
flapPointLim = 0

numberOfPoints = 0 # It is used to find the number of points in the CSV file which stores the chordline coordinates
normalizationVal = 0.01
airfoilLen = 1.0 #Length of the airfoil or the chord length
liftPerSpan = 0.0 #The output of the total lift that is generated per span is stored here
#Empties to get the camber angles
xCamber = []
yCamber = []



xCamber, yCamber = returnCamberLine(r"Airfoils\naca2412Camber.csv")
#We scale the values such that we can use it in our equations (normalization)
yCamber = [number * 1 * normalizationVal for number in yCamber]
xCamber = [number * 1 * normalizationVal for number in xCamber] 
numberOfPoints = len(xCamber)
for a in range(0, numberOfPoints):
    if xCamber[a] > (airfoilLen * (1 - percFlap/100)):
        flapPointLim = a
        break

for a in range(flapPointLim, numberOfPoints):
    yCamber[a] = yCamber[a] - (xCamber[a] - xCamber[flapPointLim])* math.tan(flapAngle)

plt.figure(figsize=(8, 6))
plt.plot(xCamber, yCamber, marker='o', linestyle='-', color='b')  # Line plot with markers

    #   Labels and title
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Plot of Y vs X')
plt.grid(True)

    # Show the plot
plt.show()
#we solve the vortex sheet equation through linear algebra A * gamma = velTerms

camberLineSlope = []
anot = 0.0
an = []
cl2 = 0

for a in  range (0, numberOfPoints-1):
    camberLineSlope.append((yCamber[a+1] - yCamber[a]) / (xCamber[a+1] - xCamber[a]))
camberLineSlope.append((yCamber[a+1] - yCamber[a]) / (xCamber[a+1] - xCamber[a]))

tetha = []

for a in range (0, numberOfPoints):
    tetha.append(np.arccos(1 - 2 * xCamber[a]))
tetha.append(tetha[a] + 0.001)

for a in range(0, numberOfPoints):
    anot = anot + camberLineSlope[a] * (tetha[a+1] - tetha[a]) 

anot = aoa - 1 / math.pi * anot

for N in range(0, numberOfPoints):
    an.append(0)
    for n in range(0, numberOfPoints):
        an[N] = an[N] + camberLineSlope[n] * math.cos(N*tetha[n]) * (tetha[n+1] - tetha[n])
    an[N] = 2 / math.pi * an[N] 
cl = math.pi * (2 * anot + an[1])
cl2 = 0
for a in range(0, numberOfPoints):
    cl2 = cl2 + camberLineSlope[a] * (math.cos(tetha[a]) - 1) * (tetha[a+1] - tetha[a])

cl2 = 2 * math.pi * (aoa + cl2 / math.pi)
print(cl)
print(cl2)