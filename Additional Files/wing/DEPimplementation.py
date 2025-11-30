"""
Assumptions-
- The airfoil is thin
- The angle of attack of the airfoil must be <<1 radian
- The jet is considered to be a rectangle that flows over the airfoil
- This rectangular jet is further approximated to be a vortex sheet

Differences from the paper implementation-
- The camber of the airfoil can be a be non symmetrical airfoil
"""
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

        
    plt.figure(figsize=(8, 6))
    plt.plot(xCamber, yCamber, marker='o', linestyle='-', color='b')  # Line plot with markers

    #   Labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Plot of Y vs X')
    plt.grid(True)

    # Show the plot
    #plt.show()
    
    return(xCamber, yCamber)

#Angle of Attack in radians
aoa = 8 * math.pi / 180

#Conditions
fsv = 20 #free stream velocity in m/s
rho = 1.225
cl = 0.0 #The coefficient of lift

#Flaps
percFlap = 30 #Percentage of the wing that is a flap
flapAngle = 40 * math.pi / 180 #angle of the 
#Jet
heightJet = 0.2 #This is responsable for setting the height of the jet (this value should include the contraction effects of the wake)
nExp = -0.1
rhoJet = 1.1
velocityJet = 50
deltaCjOverride = 0 #If this is non zero, the value of the calulated Cj will be overwritten

numberOfPoints = 0 # It is used to find the number of points in the CSV file which stores the chordline coordinates
normalizationVal = 0.01 #This is used to normalize all the x values of the airfoil such that they will be between 0 and 1
airfoilLen = 1.0 #Length of the airfoil or the chord length
liftPerSpan = 0.0 #The output of the total lift that is generated per span is stored here

#Empties to get the camber angles
xCamber = []
yCamber = []


xCamber, yCamber = returnCamberLine(r"Airfoils\NACA23015.csv") #This should be the path to where the airfoil is stored
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
tetha.append(math.pi)

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

gammaAirfoil= []
for a in range (0, numberOfPoints):
    for n in range(0, numberOfPoints):
        gammaAirfoil.append(2 * fsv * (anot * (1 + math.cos(tetha[a])) / math.sin(tetha[a]) + an[n] * math.sin(n * tetha[a])))

gammaJet = []
jetCurve = []

for a in range(0, numberOfPoints):
    jetCurve.append(yCamber[a])

mx = 1
step = 0.005
x = []
jetCurveTetha = []
Acoef = camberLineSlope[numberOfPoints-2] - aoa
yShift =  Acoef * airfoilLen / nExp - yCamber[-1]

for a in np.arange(0, mx + step, step):
    x.append( ((mx - a) / airfoilLen)  + 0.000000000000000000001)
for a in range(0, int(mx / step)):
    jetCurveTetha.append(Acoef * pow(math.e, nExp * (1 / (x[a]) - airfoilLen) / airfoilLen))
    jetCurve.append((Acoef * airfoilLen / nExp) * (pow(math.e, nExp * (1 / (x[a]) - airfoilLen) / airfoilLen)) - yShift)

totalX = []
totalX.extend(xCamber)
for a in range(0, int(mx/step)):
    totalX.append( 1 / x[a])

plt.clf()
plt.plot(totalX, jetCurve, marker='o', linestyle='-', color='b')  # Line plot with markers

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Plot of Y vs X')
plt.grid(True)

plt.show()

newGamma = []
newAn = []

deltaJ = rhoJet * velocityJet * velocityJet * heightJet - rho * fsv * fsv * heightJet
deltaCj = deltaJ / (rho * fsv * fsv * airfoilLen)
if deltaCjOverride != 0:
    deltaCj = deltaCjOverride
tethaIntegral = []
cIntegral = []
newAn.append(0)
for a in range(0, numberOfPoints):
    for xPos in range(0, int(mx/step)-1):
        if(x[xPos] > (0.0000000001 + (airfoilLen / 2 * (1 - math.cos(tetha[a]))))) or (x[xPos] < (-0.000000001 + (airfoilLen / 2 * (1 - math.cos(tetha[a]))))):
            newAn[0] = newAn[0] +  1 / math.pi * (
                deltaCj / (4 * math.pi) * (jetCurveTetha[xPos+1] - jetCurveTetha[xPos])  / (x[xPos] - airfoilLen / 2 * (1 - math.cos(tetha[a]))) * (tetha[a+1] - tetha[a])
                )
        else:
            newAn[0] = newAn[0] +  1 / math.pi * (
                deltaCj / (4 * math.pi) * (jetCurveTetha[xPos+1] - jetCurveTetha[xPos]) * (tetha[a+1] - tetha[a])
                )   
    newAn[0] = newAn[0] - 1/math.pi * camberLineSlope[a] * (tetha[a+1] - tetha[a])
newAn[0] = aoa + newAn[0]

for a in range(1, numberOfPoints):
    for b in range(0, numberOfPoints):
        for xPos in range(0, int(mx/step)-1):
            newAn.append(2 / math.pi * (
                -1 * deltaCj / ( math.pi) * (jetCurveTetha[xPos+1] - jetCurveTetha[xPos]) / (0.000000001 +  x[xPos] - airfoilLen / 2 * (1 - math.cos(tetha[b]))  ) * (tetha[b+1] - tetha[b]) * math.cos(a * tetha[b]) 
            ))
        newAn[a] = newAn[a] + 2/math.pi * camberLineSlope[b] * (tetha[b+1] - tetha[b]) * math.cos(tetha[b])

clNew = math.pi * (2 * newAn[0] + newAn[1])

print("The coefficient of lift of our blown airfoil will be =", clNew)
print("The coefficient of lift of our unblown airfoil will be =", cl)
