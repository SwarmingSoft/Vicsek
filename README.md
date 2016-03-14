# Vicsek
Vicsek simulation

tested with gcc version 4.9.2 (x86_64-win32-seh-rev2, Built by MinGW-W64 project)   
compiled with "-std=c++11", "-D_USE_MATH_DEFINES", "-DNDEBUG"  
tested with "-fexpensive-optimizations", "-O3", "-march=native"   
on Windows e.g. `g++.exe -std=c++11  -fexpensive-optimizations -O3 -march=native -DNDEBUG -IC:\voro++-0.4.6\src -IC:\eigen-eigen-10219c95fe65  -c constants.cpp -o constants.o`  
`g++.exe -std=c++11  -fexpensive-optimizations -O3 -march=native -DNDEBUG -IC:\voro++-0.4.6\src -IC:\eigen-eigen-10219c95fe65  -c data.cpp -o data.o`  
`g++.exe -std=c++11  -fexpensive-optimizations -O3 -march=native -DNDEBUG -IC:\voro++-0.4.6\src -IC:\eigen-eigen-10219c95fe65  -c tests.cpp -o tests.o`  
`g++.exe -o Vicsek.exe constants.o data.o tests.o C:\libVoro++.a`


additional dependencies:
- [voro++](http://math.lbl.gov/voro++/)
- [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)


call executable as `Vicsek.exe data output_mode tskip tevery tconsecutive mode tend N L dt v0 J eta seed geometry geometry_param`

where  
output (data): output file "data.txt"  
output_mode (0): 0 = data.txt and data.extras.txt with extra information, 1: only data.txt, 2: only data.extras.txt  
tskip (20000): skip the forst tskip timesteps  
tevery (4000): output only every tevery timesteps  
tconsecutive (1): 0: no additional timesteps, 1: for each output also output the following timestep  
mode (vicsek): vicsek: normal Vicsek, linearvicsek: linear Vicsek, linearvicseknoisetransform: transforms the noise into the eigen space, subtract random walk component and transforms back (untested)  
tend(420002): timesteps after which the ismulation terminates  
N (256): number of particles  
L (16): length of the square simulation box  
dt (0.01): euler integrator timestep  
v0 (1.): constant particle absolute velocity  
J (1.): coupling constant  
eta (0.2): noise amplitude (uniform distribution)  
seed (726): random number generator seed  
geometry (metric): metric: normal metric topology, topological: symmetrised nearest neighbour topology, voronoi: Voronoi topology  
geometry_param (0.5): radius in metric case, number of neighbours in topological case, leave this parameter out for voronoi