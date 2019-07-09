#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 15:52:26 2019

@author: whitefar
"""
#
from shapely.geometry import Point, LineString, MultiPoint
import numpy as np
import geopandas as gpd
import os
import pygmsh.built_in as pyg
import pygmsh
import meshio
import matplotlib.pyplot as plt
#


class mesh:       
        
    def load_perimeter(self,perimeter_file):
        """
        Input filetype must be readable by geopandas .read_file, and must be a linestring, or polygon. E.g. .shp or .gpx
        
        see import fiona; help(fiona.open) for info.
        """
        self.perimeter_file = perimeter_file
            
        #Return a GeoDataFrame object using geopandas
        self.perimeter_gdf = gpd.read_file(perimeter_file)
        
        if self.perimeter_gdf.iloc[0].geometry.type =='Polygon':
            self.geometry = self.perimeter_gdf.iloc[0].geometry
            self.perimeter = LineString(list(self.geometry.exterior.coords))
            self.perimeter_coords = list(self.geometry.exterior.coords)
            self.perimeter_array = np.array(list(self.geometry.exterior.coords))
            self.num_points =  len(list(self.geometry.exterior.coords))
        elif self.perimeter_gdf.iloc[0].geometry.type =='LineString':
            self.geometry = self.perimeter_gdf.iloc[0].geometry
            self.perimeter = self.geometry
            self.perimeter_coords = list(self.perimeter.coords)
            self.perimeter_array = np.array(list(self.perimeter.coords))
            self.num_points = len(list(self.perimeter.coords))
        elif self.perimeter_gdf.iloc[0].geometry.type == 'Point':
            self.perimeter = LineString([self.perimeter_gdf.iloc[i].geometry.coords[:][0] for i in range(len(self.perimeter_gdf))])
            self.perimeter_coords = list(self.perimeter.coords)
            self.perimeter_array = np.array(list(self.perimeter.coords))
            self.num_points = len(list(self.perimeter.coords))
        else:
            raise ValueError("Input must have Point, LineString or Polygon geometry")
        
        print(self.perimeter_gdf.head())
        self.perimeter_gdf.plot()
        
        
            
    
    def set_resolutions(self,resolutions):
        """
        specify the output mesh grid spacing. Should be a single interger/float 
        or a list. Units are the same as perimeter coordinates.
        """
        
        if type(resolutions)=='int' or type(resolutions)=='float':
            self.resolutions = [resolutions]
        else: 
            self.resolutions = resolutions
            
                
        
    
    def make_perimeters_newres(self):
        """
        must have run set_resolutions().
        make_perimeters_newres makes a list of perimeter files each with resolution/grid spacing
        as specified by set_resolutions()
        """

        
        # this writes
        
        
        
        self.perimeters_newres = []
        
        
        for distance in self.resolutions:
            
            #write the points using shapely
            list_points = []
            ## set the current distance to place the point
            current_dist = distance
            ## make shapely MultiLineString object
            shapely_line = self.perimeter
            ## get the total length of the line
            line_length = shapely_line.length
            ## append the starting coordinate to the list
            
            list_points.append(Point(list(shapely_line.coords)[0]))
            ## https://nathanw.net/2012/08/05/generating-chainage-distance-nodes-in-qgis/
            ## while the current cumulative distance is less than the total length of the line - last step
            while current_dist < (line_length - distance):
                ## use interpolate and increase the current distance
                list_points.append(shapely_line.interpolate(current_dist))
                current_dist += distance
           
            self.perimeters_newres.append(LineString([(lp.x,lp.y) for lp in list_points]))
            
            #plt.plot(self.perimeters_newres[-1].xy[0],self.perimeters_newres[-1].xy[1],'o')
        
        self.perimeters_newres_coords = [list(self.perimeters_newres[i].coords) for i, _ in enumerate(self.perimeters_newres) ]
        self.perimeters_newres_arrays = [np.array(list(self.perimeters_newres[i].coords)) for i, _ in enumerate(self.perimeters_newres) ]
            
    def plot(self,resolution_index=None):
        """
        plots the original perimeter with points from a particular resolution.
        Specify the index of the resolution in the .resolution list.
        """

        plt.plot(self.perimeter.xy[0],self.perimeter.xy[1])
        if resolution_index != None:
            plt.plot(self.perimeters_newres[resolution_index].xy[0],self.perimeters_newres[resolution_index].xy[1],'o')
            print("Plotted with points {} apart".format(self.resolutions[resolution_index]))
    
    
    def make_meshes(self):
        """ Makes a 2D mesh from a border which ElmerGrid can interpret
        
        INPUT: a border .xy file, distanc is characteristic distance between points
        
        OUTPUT: a .geo file of the border (gmsh geo object format)
           and gmsh mesh .msh file for use in ElmerGrid
        """
        
        mesh_file_type='msh'
        input_root, input_ext = os.path.splitext(self.perimeter_file)
        self.output_filenames = []

        #for each different resolution
        for i, distance in enumerate(self.resolutions):
            input_array = np.array(self.perimeters_newres[i])
            length, q = input_array.shape
            input_array = np.hstack((input_array, np.zeros([length,1]) ))
            
            #load np array to pygmsh geometry class
            geom = pyg.Geometry()
            poly = geom.add_polygon(input_array,distance)
            #poly = geom.add_polygon(input_array,dist)
            #geom.add_physical_surface(poly,'poly')
            
            #make a mesh out of it

            output_filename = input_root +"_"+ str(distance)
            mesh = pygmsh.generate_mesh(geom,verbose=False,dim=2,geo_filename=output_filename+".geo",mesh_file_type=mesh_file_type)
            self.output_filenames.append(output_filename)
            
            # write the mesh to gmsh msh format using meshio
            meshio.write_points_cells(
                output_filename+"."+mesh_file_type,
                mesh.points,
                mesh.cells,
                point_data=mesh.point_data,
                cell_data=mesh.cell_data,
                field_data=mesh.field_data,
                file_format='gmsh2-ascii'
                )
            
            # write the mesh to vtu for viewing in paraview
            meshio.write_points_cells(
                output_filename+".vtu",
                mesh.points,
                mesh.cells,
                point_data=mesh.point_data,
                cell_data=mesh.cell_data,
                field_data=mesh.field_data
                )
   
            print(input_root+" written as "+output_filename+".geo and "+output_filename+".msh and "+output_filename+".vtu")
        return 


    def export(self,filetype=".shp",driver='ESRI Shapefile'):
        """
        indicies = index of mesh res can be one, or a list or "all"
        """
        
        
        input_root, input_ext = os.path.splitext(self.perimeter_file)
        out_gdf = self.perimeter_gdf
                
        for index, resolution in enumerate(self.resolutions):
            out_gdf.loc[0, 'geometry'] = self.perimeters_newres[index]
            out_gdf.to_file(input_root+"_"+str(resolution)+filetype,driver='ESRI Shapefile')
            print("exported "+input_root+"_"+str(resolution)+filetype+" to {}".format(os.getcwd()))
            
            
    def print_elmergrid_command(self,run_as="print"):
        """
        """
        
        if run_as=="script":
                
            with open('run_elmergrid.sh','w') as runfile:
                runfile.write('#!/bin/bash\n')
                for name in self.output_filenames:
                    runfile.write("ElmerGrid 14 2 {} -autoclean\n".format(name))
            print("the bash script 'run_elmergrid.sh' has been written to {}. This script processes all mesh files with ElmerGrid.".format(os.getcwd()))
        else:
            print(" ")
            print(" ")
            print(" ")
            print(" ")
            print("Printing a command to run in bash shell. This command processes all mesh files with ElmerGrid.".format(os.getcwd()))
            print(" ")
            print("The command is below:")
            print(" ")
            print(" ")
            print(" ")
            print(" ")
            print("cd {}".format(os.getcwd()))
            for name in self.output_filenames:
                print("ElmerGrid 14 2 {} -autoclean".format(name))
            print(" ")
            print(" ")
            print(" ")
            print(" ")

                
    def run_all(self):
        self.spaced_points()
        self.make_meshes()
        self.export()
        self.print_elmergrid_script()
        
    
    
    def find_corners(self,indicies):
        
        plt.plot(self.perimeter.xy[0],self.perimeter.xy[1])
        
        if len(indicies)==4:
            plt.plot(self.perimeter.xy[0][indicies[0]],self.perimeter.xy[1][indicies[0]],'ro')
            plt.plot(self.perimeter.xy[0][indicies[1]],self.perimeter.xy[1][indicies[1]],'go')
            plt.plot(self.perimeter.xy[0][indicies[2]],self.perimeter.xy[1][indicies[2]],'o',color='orange')
            plt.plot(self.perimeter.xy[0][indicies[3]],self.perimeter.xy[1][indicies[3]],'bo')
            plt.grid()
            print("Perimeter has {} points.".format(self.num_points))
            print("Red dot is at {}. Green at {}. Orange at {}. Blue at {}.".format(self.perimeter.coords[indicies[0]]\
                  ,self.perimeter.coords[indicies[1]],self.perimeter.coords[indicies[2]],self.perimeter.coords[indicies[3]]))
                    
            self.found_corners = [self.perimeter.coords[indicies[0]]\
                  ,self.perimeter.coords[indicies[1]],self.perimeter.coords[indicies[2]],self.perimeter.coords[indicies[3]]]
        else:
            cmap=plt.get_cmap('jet',self.num_points)
            self.found_corners = []
            for index in indicies:
                plt.plot(self.perimeter.xy[0][index],self.perimeter.xy[1][index],'o',color=cmap(index/self.num_points),label=str(index))
                self.found_corners.append(self.perimeter.coords[index])
            plt.legend()
            print("Perimeter has {} points.".format(self.num_points))
            print("Each index provided corresponds to a coloured point, shown in the legend")
            
        self.found_corner_indicies = indicies
        
    
    def set_corners(self,corners=[]):
        """
        will take self.find_corners
        alternatively, input your own coordinates.
        """
        if corners == []:
            corners=self.found_corners
            self.corner_indicies = self.found_corner_indicies
        else:
            self.corner_indicies = []
            
            for corner in corners:
                corner_index = np.argmin(np.sum((self.perimeter_array- corner)**2, axis=1))
                self.corner_indicies.append(corner_index) 
        
        plt.plot(self.perimeter.xy[0],self.perimeter.xy[1])
        for corner in corners:
            plt.plot(corner[0],corner[1],'ro')
        self.corners = corners
        
        
        
        print("corners that divide boundaries are set as the red dots.")
        

    
    def label_boundaries(self, labels = [3,1,3,2]):
        """label elmer boundaries for use in ElmerSolver
         
            
            
        """
        self.boundary_labels = labels
        self.corners_newres_indicies = []
        
        #for each different resolution
        for i,mesh_name in enumerate(self.output_filenames):
            
            
        #DOUBLE CHECK THE DIRECTION OF THE PATH
        
        #mesh_folder,perimeter_xy ,distance,
            distance = self.resolutions[i]
        
            mesh_file = mesh_name+"/mesh.boundary"
            
            boundary_array = np.loadtxt(mesh_file)
            (nrows,ncols) = boundary_array.shape
            
            #FIRST
            # Get the resolution specific indicies for the corners - verticies where the different edge labels are divided
            corners = np.array(self.corners)
            new_corner_indicies = []
            
            for corner in corners:
                new_corner_index = np.argmin(np.sum((self.perimeters_newres_arrays[i]- corner)**2, axis=1))
                new_corner_indicies.append(new_corner_index)  
                
            self.corners_newres_indicies.append(new_corner_indicies)
            
           #SECOND make a vector of boundary labels to match the corner indicies
           # First four sides [start-corner1, corner1-corner2,corner2-corner3,corner3-corner4]
            boundary_labels = []
            corner_j_prev = 0 
           
            for j, corner_j in enumerate(new_corner_indicies):
                len_side = int(corner_j - corner_j_prev)
                if len_side<0:
                    len_side = len(self.perimeters_newres_coords[i]) + len_side
                add_labels = np.ones((len_side),dtype=int)*labels[j]
                boundary_labels = np.append(boundary_labels,add_labels)
                corner_j_prev = corner_j
           
            #Fifth side to return to start [corner4 - start]
            #Shouldnt run if perimeter starts on a corner
            if nrows -  boundary_labels.shape[0] != 0:
                add_labels = np.ones((nrows -  boundary_labels.shape[0]),dtype=int)*labels[0]
                boundary_labels = np.append(boundary_labels,add_labels)
            
            # substitute these labels in the elmer boundary file 
            boundary_array[:,1] = boundary_labels
            np.savetxt(mesh_file,boundary_array,fmt='%i')
              
            print("mesh.boundary overwritten with new boundary labels for grid spacing of {}".format(distance))
            
        plt.plot(self.perimeter.xy[0],self.perimeter.xy[1])
        plt.title("Boundary labels")
        for corner in self.corners:
            plt.plot(corner[0],corner[1],'ro')
        for i, label in enumerate(labels):
            half_index = int(self.corner_indicies[i]/2+self.corner_indicies[i-1]/2)
            if half_index>self.corner_indicies[i]:
                half_index = int(self.corner_indicies[i]/2+(self.corner_indicies[i-1]-self.num_points)/2)
            plt.text( self.perimeter_coords[half_index][0],\
                      self.perimeter_coords[half_index][1],str(label),fontsize=15,weight='bold' )
    
    
    def plot_boundaries(self, index):
        
        new_res_perimeter = self.perimeters_newres_arrays[index]
        corner_indicies = self.corners_newres_indicies[index]
        
        plt.plot(self.perimeter_array[:,0],self.perimeter_array[:,1])
        plt.plot(new_res_perimeter[:,0],new_res_perimeter[:,1],'bo')
        plt.title("Boundary labels")
        for i in corner_indicies:
            plt.plot(new_res_perimeter[i,0],new_res_perimeter[i,1],'ro',markersize=15)
        for i, label in enumerate(self.boundary_labels):
            half_index = int(self.corner_indicies[i]/2+self.corner_indicies[i-1]/2)
            if half_index>self.corner_indicies[i]:
                half_index = int(self.corner_indicies[i]/2+(self.corner_indicies[i-1]-self.num_points)/2)
            plt.text( self.perimeter_coords[half_index][0],\
                      self.perimeter_coords[half_index][1],str(label),fontsize=20,weight='bold',color='g' )
        print("Red dots separate green labeled boundaries")





# perimeter_file = 'example_perimeter/brewster.shp'
# msh = mesh()
# msh.load_perimeter(perimeter_file)
# msh.set_resolutions([100,200])
# msh.make_perimeters_newres()
# msh.plot(1)
# #msh.perimeters_newres[2]
# msh.make_meshes()

# msh.print_elmergrid_command()
# #msh.run_all()
# msh.find_corners([100,200,300,400,500])
# msh.set_corners()
# msh.label_boundaries([1,2,3,4,5])
        
