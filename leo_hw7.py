import numpy as np
import matplotlib.pyplot as plt

n = np.linspace(0, 7.9, 200)
v = np.zeros_like(n)

for i, val in enumerate(n):
    v[i] = np.degrees(np.atan(2 * np.tan(np.radians(val))))

plt.plot(n, v)
plt.xlabel("Look Angle (η) (degrees)")
plt.ylabel("Vertex Angle of Cone (degrees)")
plt.show()