# PostProcessing Ligghts

This software provides many post-processing functions to handle data simulated using the Discrete Element Modelling software Liggghts. It started quite a while ago and has accumulated years of code, so no everything is relevent of working. Nevertheless, it can be quite powerful and useful.

# Getting started
## Compiling
These instruction are provided with a linux os in mind. Compiling on other operating systems is possible, but your milage may vary. 

- Clone the repository: `git clone https://github.com/Franzzzzzzzz/PostProcessingLiggghts`

- Move there: `cd PostProcessingLiggghts`

- Create a build folder: `mkdir build && cd build`

- Prepare the compilation: `cmake ..`

  - MATLAB support: if you want to support writing the data in Matlab format, cmake needs to be able to find the matlab library. This is done by adjusting the variable `Matlab_ROOT_DIR` to your system installation. 
  
  - Other useful library: libtif, zlib, boost. Some of these may allow compilation even when not found, other won't.
  
  - Language: for historical reason, the default language of the software is a mix of mostly French and English. If you want English only, set the `language` variable in cmake to `en_AU.UTF-8`.
  
- Build: `make`
- You should now have an executable called `PostProcessing` in the build folder. 



## Using PostProcessing
The program runs in command line only. The last argument is always a dump from Liggghts. The argument before the last one can be a dump from ligghts of the particle positions if the last argument is a dump of contact forces. 

Examples: 

- Get help: `PostProcessing --help null`

- Create vtk files from dump files: `PostProcessing --dump2vtk --W '[vtk]def,v,type,radius' dump.test`

- Create a vtk with the contact forces: `PostProcessing --dump2vtk --chainforce --W '[vtk]def,f' dump.test dump.force`

- Create a vtk with orientation information for superquadric/ellipsoidal particles, to use with the TensorGlyph filter in Paraview (requires to export the quaternion components quat1-4): `PostProcessingGamma --dump2vtk --superquadric 0.001 0.001 0.002 --W '[+vtk]def,v,orient,radius' dump.test`

- Coarse-grain particle data: `PostProcessing --coarse-graining 3 4 5 --W "[vtk]v,phi" dump.test`

- Coarse-grain contact data: `PostProcessing --coarse-graining 3 4 5 --use-box -0.2 0.2 -0.1 0.1 1.2 1.5 --chainforce --W "[vtk]sigmatot" --w/kinetic dump.particles dump.contact`

