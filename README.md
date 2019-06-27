# elmer_mesh_maker

This module makes 2D meshes of different resolutions for ElmerIce, with customised labeled boundaries. ElmerIce will extrude the mesh into 3D.

See jupyter notebook tutorial.

Process summary:
1. msh = mesh()
2. msh.load_perimeter(perimeter_file)
3. msh.set_resolutions(list_of_resolutions)
4. msh.make_perimeters_newres()
5. msh.make_meshes()
6. (optional) msh.export()  
7. msh.print_elmergrid_command()
8. msh.set_corners(list_of_coordinates)
9. msh.label_boundaries(list_of_labels)
        
jupyter notebook tutorial:

# # Create a mesh with labeled boundaries for ElmerIce
# 
# This module makes 2D meshes of different resolutions for ElmerIce, with customised labeled boundaries.  ElmerIce will extrude the mesh into 3D.
# 
# 
# ### Required python modules:
# 
# * meshio
# * pygmsh
# * shapely
# * geopandas 
# * fiona 
# * matplotlib
# * numpy
# * os
# 
# 
# ## Process:
# 
# **Making a mesh**
# 1. Load mesh maker package
# 2. Load glacier perimeter file
# 3. Make 1D perimeters with different resolutions (grid spacing)
# 4. (optional) Export the new resolution perimeters.
# 5. Make 2D meshes for all set resolutions
# 6. Use ElmerGrid to make the meshes suitable for Elmer
# 
# **Labeling the boundaries**
# 
# 1. Set the glacier corners
# 2. Label the glacier boundaries for ElmerIce to read
# 
# ------------------------------------------------------------------
# 
# 
# # Making a mesh
# 
# ## 1. Load mesh maker package

# In[1]:


import os
import sys

#location of elmer_mesh_maker package
package_location = 'elmer_mesh_maker'
sys.path.append(os.path.abspath(package_location))

import elmer_mesh_maker as el


# ## 2. Load input glacier perimeter

# In[2]:


#create mesh object
msh = el.mesh() 

#change directory to a folder to read perimeter file from and write meshes to
folder = 'elmer_mesh_maker' 
os.chdir(folder)

perimeter_file = 'brewster.shp'

msh.load_perimeter(perimeter_file)


# ## 3. Make perimeters with different resolutions (grid spacing)
# The mesh maker gmsh works by taking a 1D perimeter and extrapolating a 2D mesh. The grid spacing in the perimeter will determine the grid spacing of the mesh.
# 
# First set resolutions, then make the perimeters

# In[3]:


msh.set_resolutions([50,100,200,300]) #units are the same as input perimeter coordinates

msh.make_perimeters_newres()

#put the index of the resolution you wish to plot
msh.plot(0)


# In[4]:


msh.plot(3)


# ## 4. (optional) Export the new resolution perimeters.

# In[5]:


msh.export()


# ## 5. Make meshes for all set resolutions
# 
# make_meshes() uses pygmsh and meshio which use gmsh to make the meshes in a gmsh .msh format. 

# In[6]:


msh.make_meshes()


# ## 6. Use ElmerGrid to make the meshes suitable for Elmer
# 
# Use the following cell to make a command you can copy and paste into the bash shell to run ElmerGrid over the meshes. 
# 
# ElmerGrid to modifys the mesh to an Elmer format. More info [here](http://www.elmerfem.org/blog/preprocess/mesh-generation-software-used-with-elmer/)
# 
# Alternatively put ```msh.print_elmergrid_command('script') ``` to print a bash script to file which will run the same command.

# In[7]:


msh.print_elmergrid_command()        


# List all the mesh names. The Elmer-ready meshes show up as a folder with mesh.boundary mesh.header mesh.node and mesh.elements files.

# In[8]:


msh.output_filenames


# # Labeling the mesh boundaries
# 
# To use ElmerIce, the mesh must have labeled boundaries. ElmerGrid automatically labels the boundaries, but you may want them in separate places.
# 
# ## 1. Set coordinates for the corners
# 'Corners' are points which separate the desired boundaries.

# In[11]:


msh.set_corners([(1314274.9987104018, 5112882.022700516),
 (1314266.9226336845, 5113679.104938047),
 (1315070.9465903954, 5114035.080386854),
 (1314553.726167038, 5112929.071662091)])


# ### OR
# If you cant be arsed opening a map to find the coordinates, use the following command to find the corner coordinates by guessing and checking the corner point indicies

# In[12]:


msh.find_corners([4,46,614,680])
msh.found_corners
msh.set_corners()


# ## 2. Label the glacier boundaries for ElmerIce to read
# Input a list of integers

# In[13]:


msh.label_boundaries([2,3,1,3])


# Here, I wanted the walls to both be 3, the glacier head one, and the terminus 2.
# 
# Plot the boundaries to check that corners and labels are in the right place for different resolutions.

# In[14]:


msh.plot_boundaries(0)


# In[15]:


msh.plot_boundaries(3)


# # END
# The mesh is ready for elmer
# 
# To visualise the mesh, open ElmerGUI and load the mesh folder. Select View -> Numbering -> Boundary Index to check the boundary labeling is right

# ### Other useful commands:

# In[16]:


msh.perimeter #the shapely linestring of input perimeter


# In[17]:


msh.perimeters_newres  # a list of the shapely linestrings of set resolution perimeters eg:
msh.perimeters_newres[3]


# In[18]:


msh.found_corners #the coordinates of corners found with find_corners()

