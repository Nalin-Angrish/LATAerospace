'''
Paper implementation of the DEP with the following considerations-
- The airfoil is thin
- It is a symmetric flat airfoil with a flap at the end which is responsible for the deflection of the jet of aiir
- The angle of attacks must be low (Turbulent air has not been modelled
- A rectangular Jet has been assumed to flow all over the surface of the airfoil
- The contraction of this jet has not yet been taken into account
'''

import math

velJet = 50
heightJet = 0.2
rhoJet = 2
rhoinf = 1
chord = 1
Vinf = 20
aoa = 10 * math.pi / 180
jetDeflection = 5 * math.pi / 180

jDash = rhoJet * velJet * velJet * heightJet
cj = jDash / (rhoinf * Vinf * Vinf * chord)
cl = 2 * math.pi * (1 + 0.151 * math.sqrt(cj) + 0.219 * cj) * aoa + 2 * math.sqrt(math.pi * cj) * math.sqrt(1 + 0.151 * math.sqrt(cj) + 0.139 * cj) * jetDeflection


print(cl)
