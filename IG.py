# -*- coding: utf-8 -*-
"""
Created on Fri May 14 07:39:02 2021

@author: Magda
"""

import geopandas as gpd
import numpy as np
import shapely
import matplotlib.pyplot as plt

gdf = gpd.read_file("PD_STAT_GRID_CELL_2011.shp")
woj = gpd.read_file("Wojewodztwa")
gdf.to_crs("EPSG:4326")
gdf['centroid'] = gdf.centroid



gdf.to_crs("EPSG:4326")
xmin, ymin, xmax, ymax= [13 ,48 , 25, 56]
n_cells=30
cell_size = (xmax-xmin)/n_cells
grid_cells = []
for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
    for y0 in np.arange(ymin, ymax+cell_size, cell_size):

        x1 = x0 - cell_size
        y1 = y0 + cell_size
        grid_cells.append( shapely.geometry.box(x0, y0, x1, y1) )
cell = gpd.GeoDataFrame(grid_cells, columns=['geometry'])
ax = gdf.plot(markersize =.1, figsize =12, column ='TOT', cmap ='jet')
plt.autoscale(False)
cell.plot(ax = ax, facecolor = "none", edgecolor = 'grey')
ax.axis("off")
merged = gpd.sjoin(gdf, cell, how='left', op='within')
dissolve = merged.dissolve(by="index_right", aggfunc="sum")
cell.loc[dissolve.index, 'TOT_0_14'] = dissolve.TOT_0_14.values
ax = cell.plot(column='TOT_0_14', figsize=12, cmap='viridis', vmax=700000, edgecolor="grey", legend = True)
ax.set_axis_off()
plt.title('liczba ludności z przedzialu 0-14 w siatce')
merged = gpd.sjoin(gdf, cell, how='left', op='within')
dissolve = merged.dissolve(by="index_right", aggfunc ="sum")
cell.loc[dissolve.index, 'TOT_15_64'] = dissolve.TOT_0_14.values
ax = cell.plot(column='TOT_15_64', figsize = 12, cmap ='viridis', vmax = 700000, edgecolor = "grey", legend = True)
ax.set_axis_off()
plt.title('liczba ludności z przedzialu 15-64 w siatce')
merged = gpd.sjoin(gdf, cell, how='left', op ='within')
dissolve = merged.dissolve(by ="index_right", aggfunc ="sum")
cell.loc[dissolve.index, 'TOT_65__'] = dissolve.TOT_65__.values
ax = cell.plot(column ='TOT_65__', figsize =12, cmap='viridis', vmax = 700000, edgecolor="grey", legend = True)
ax.set_axis_off()
plt.title('liczba ludności z przedzialu >65 w siatce')
merged = gpd.sjoin(gdf, cell, how ='left', op ='within')
dissolve = merged.dissolve(by ="index_right", aggfunc ="sum")
cell.loc[dissolve.index, 'TOT_65__'] = dissolve.MALE_0_14.values + dissolve.MALE_15_64.values + dissolve.MALE_65__.values
ax = cell.plot(column='TOT_65__', figsize = 12, cmap = 'viridis', vmax=700000, edgecolor="grey", legend = True)
ax.set_axis_off()
plt.title('liczba ludności z przedzialu >65 w siatce')
merged = gpd.sjoin(gdf, cell, how ='left', op = 'within')
dissolve = merged.dissolve(by ="index_right", aggfunc ="sum")
cell.loc[dissolve.index, 'TOT_65__'] = dissolve.FEM_0_14.values + dissolve.FEM_15_64.values + dissolve.FEM_65__.values
ax = cell.plot(column='TOT_65__', figsize = 12, cmap = 'viridis', vmax = 700000, edgecolor = "grey", legend = True)
ax.set_axis_off()
plt.title('liczba ludności z przedzialu >65 w siatce')

ax = cell.plot(column ='TOT', figsize = 12, cmap ='viridis', vmax = 700000, edgecolor ="grey", legend = True)
ax.set_axis_off()
plt.axis('equal')
plt.title('Liczba ludności w siatce')
merged = gpd.sjoin(gdf, cell, how='left', op='within')
dissolve = merged.dissolve(by ="index_right", aggfunc ="sum")
cell.loc[dissolve.index, 'TOT'] = dissolve.TOT.values
 #jnrke
