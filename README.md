# connectivityMap
Calculate and create connectivity map and matrix from lagrangian trajectory model
This project is based on Opendrift, but can be used with different models
To create a connectivity map and matrix run the following functions:
1. `run_for_connectivity_polygons` \ `run_for_connectivity_domain` - to run the lagrangian model and create the sectors, this will return 
    1. `first_lonlats` - `{'lon':lons,'lat':lats}`, the origin of the elements of the model
    2. `last_lonlats` - `{'lon':lons, 'lat':lats}`, the last location of the elements of the model
    3. `supply` - the number of elements supplied by each sector
    4. `sectors` - a list containing (lon, lat) points of the verticies of the sectors
  This can be step skiped if you already have this parameters
2. `calculate_connectivity_matrix(first_lonlats, last_lonlats, supply, sectors)` - to create the connectivity matrix from the lagrangian model. This will return `C` - the connectivity matrix
3. `categorize_sectors(sectors, C)` - to categorize the sectors by connectivity categories
4. `plot_connectivity_pylygons(sectors, sector_categories)` \ `plot_connectivity_domain(sectors, sector_categories)` - this will plot the `sectors` and color each sector with an appropriate color to its value from `sector_categories` (i.e. `sector[i]` will be plotted with the color for the number at `sector_categories[i]`)

![connectivity matrix](https://github.com/yoelbassin/connectivityMap/blob/master/data/pictures/connectivity_matrix.png?raw=true)


