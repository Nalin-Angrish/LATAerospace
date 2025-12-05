import math
import numpy as np
import matplotlib.pyplot as plt

chord_A = 2
chord_alpha = 1
chord_beta = 9
chord_fun = lambda r: chord_A * (r**chord_alpha) * ((1 - r)**chord_beta)

radii = np.linspace(0.05, 0.4, 100)
chord = [chord_fun(r) for r in radii]
plt.plot(radii, chord)
plt.xlabel("Radius (m)")
plt.ylabel("Chord (m)")
plt.title("Blade Chord Distribution")
plt.grid(True)
plt.show()