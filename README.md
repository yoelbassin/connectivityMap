# connectivityMap
Calculate and create connectivity map and matrix from lagrangian trajectory model
This project is based on Opendrift, but can be used with different models
To create a connectivity map and matrix run the following functions:
1. `Sector.create_sectors_domain` \ `Sector.create_sectors_domain` - To create the sectors for the model to analize
2. `model = Connectivity(lons, lats, sectors)` using the lons and lats of the lagrangian trajectories, and sectors generated from (1).
3. `model.plot(model.n_settle_values, title=Nsettle)` will plot the model with its 'n_settle' values.
<p align="center">
    <img src="https://github.com/yoelbassin/connectivityMap/blob/master/data/pictures/n_settle.png" width=400>
</p>

