import numpy as np
import ase.atoms
import ase.io
import ase.visualize

cell = np.array([[0, 0, 0], [1, 1, 0], [1, 0, 1], [0, 1, 1]])*.5
a = 4.45
N = 6

images = np.array([(i % N, j % N, k % N) for i in range(N)
                   for j in range(N) for k in range(N)])

x = np.zeros((4*N**3, 3))

for i, img in enumerate(images):
    x[i*4:(i+1)*4, :] = a*(cell+img)

atoms = ase.Atoms(symbols=f'Ne{len(x)}', positions=x,
                  # velocities=np.random.normal(0, 3.51e-4*1**0.5, (len(x), 3)),
                  cell=(N*a*np.identity(3)).round(2), pbc=True)
ase.io.write('positions.in', atoms, format='extxyz')
# ase.visualize.view(atoms)
