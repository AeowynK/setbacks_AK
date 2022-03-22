import geopandas as gpd
import matplotlib.pyplot as plt
hood = gpd.read_file('.\data\shapefiles\\neighborhood.shp')
hood.plot(alpha=0.5)
plt.show()

#boundary = gpd.read_file('.\data\shapefiles\\boundary.shp')
#boundary.boundary.plot(alpha=0.5)
#plt.show()