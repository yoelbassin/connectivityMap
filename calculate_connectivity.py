# This file is part of ConnectivityMap module
#
# ConnectivityMap is a free module for creating a connectivity map 
# and matrix from a lagrangian trajectory model
#
# Copyright 2021, Erick Fredj and Yoel Bassin, JCT, Israel


import numpy as np
from datetime import timedelta
import pyproj
import copy
import cartopy
import matplotlib.pyplot as plt
import matplotlib
from opendrift.readers import reader_global_landmask
from shapely.geometry.polygon import Polygon
import cartopy.crs as ccrs

"""
This file contains the functions of the ConnectivityMap module.
"""


# Preprocessing and helpful functions

def is_point_in_area(lon_list, lat_list, point):
    """
    checks if a point is inside of a given polygon

    :param lon_list: list of longitudes of the  polygon
    :param lat_list: list of latitudes of the polygon
    :param point: coordinate (lon, lat) of the point

    :return: boolean
    """
    if type(point) != tuple or len(point) != 2:
        raise ValueError('A vertex is a tuple of length 2')
    nvert = len(lon_list)
    vertx = lon_list
    verty = lat_list
    testx, testy = point[0], point[1]
    c = 0
    j = nvert - 1
    for i in range(0, nvert):
        if (((verty[i] > testy) != (verty[j] > testy)) and
                (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i])):
            c = not (c)
        j = i
    return c


def sectors_from_shapefile(shapefile):
    """
    creates a list of polygon sectors from a shapefile
    from opendrift

    :param shapefile: shapefile path

    :return: list of polygon sectors [(lon_list, lat_list), (lon_list, lat_list), ...]
    """
    sectors = []
    layername = None
    featurenum = None
    number = 100
    try:
        from osgeo import ogr, osr
    except Exception as e:
        raise ValueError('OGR library is needed to read shapefiles.')

    targetSRS = osr.SpatialReference()
    targetSRS.ImportFromEPSG(4326)
    try:
        s = ogr.Open(shapefile)
    except:
        s = shapefile

    for layer in s:
        if layername is not None and layer.GetName() != layername:
            print('Skipping layer: ' + layer.GetName())
            continue
        else:
            print('Seeding for layer: %s (%s features)' %
                  (layer.GetDescription(), layer.GetFeatureCount()))

        coordTrans = osr.CoordinateTransformation(layer.GetSpatialRef(),
                                                  targetSRS)

        if featurenum is None:
            featurenum = range(1, layer.GetFeatureCount() + 1)
        else:
            featurenum = np.atleast_1d(featurenum)
        if max(featurenum) > layer.GetFeatureCount():
            raise ValueError('Only %s features in layer.' %
                             layer.GetFeatureCount())

        # Loop first through all features to determine total area
        layer.ResetReading()
        area_srs = osr.SpatialReference()
        area_srs.ImportFromEPSG(3857)
        areaTransform = osr.CoordinateTransformation(
            layer.GetSpatialRef(), area_srs)

        areas = np.zeros(len(featurenum))
        for i, f in enumerate(featurenum):
            feature = layer.GetFeature(f - 1)  # Note 1-indexing, not 0
            if feature is not None:
                gom = feature.GetGeometryRef().Clone()
                gom.Transform(areaTransform)
                areas[i] = gom.GetArea()
        total_area = np.sum(areas)
        layer.ResetReading()  # Rewind to first layer
        print('Total area of all polygons: %s m2' % total_area)
        # Find number of points per polygon
        numbers = np.round(number * areas / total_area).astype(int)
        numbers[numbers.argmax()] += int(number - sum(numbers))

        for i, f in enumerate(featurenum):
            feature = layer.GetFeature(f - 1)
            if feature is None:
                continue
            geom = feature.GetGeometryRef()
            try:
                geom.Transform(coordTrans)
            except Exception as e:
                print('Could not transform coordinates:')
                print(e)
                pass
            # b = geom.GetBoundary()
            # if b is not None:
            #    points = b.GetPoints()
            #    lons = [p[0] for p in points]
            #    lats = [p[1] for p in points]
            # else:
            # Alternative if OGR is not built with GEOS support
            r = geom.GetGeometryRef(0)
            lons = [r.GetY(j) for j in range(r.GetPointCount())]
            lats = [r.GetX(j) for j in range(r.GetPointCount())]
            if not lons or not lats:
                continue
            sectors.append((lons, lats))
    return sectors


# Seeding elements and running simulation

def run_for_connectivity_polygons(opendrift_model, time=None, steps=None, duration=None,
                                  end_time=None, time_step=None, depth=0, sectors=[], supply_ammount=10,
                                  time_step_output=None):
    """
    seed and run OpenDrift model before creation of a connectivity map using polygon sectors

    :param opendrift_model: opendrift model for the simulation
    :param time: datenum or list.
        one element - The time at which particles are seeded/released.
        list with two elements - elements are seeded continuously from start/first to end/last time.
        list with more than two elements - the number of elements is equal to len(time) and are seeded as a time series.
    :param steps: integer, maximum number of steps. End of simulation will be self.start_time + steps*self.time_step
    :param duration: timedelta defining the length of the simulation
    :param end_time: datetime object defining the end of the simulation
    :param time_step: interval between particles updates, in seconds or as timedelta. Default: 3600 seconds (1 hour)
    :param depth: depth at which the elements are seeded
    :param d_supply: the amount of elements to supply at each sector
    :param time_step_output: Time step at which element properties are stored and eventually written to file.
    Timedelta object or seconds. Default: same as time_step, meaning that all steps are stored

    :return: {'lons': initial_lons, 'lats': initial_lats}, {'lons': final_lons, 'lats': final_lats},
    supply, sectors - the origin location of elements, the amount supplied by each sector and the
    coordinates of the sectors
    """

    ########################
    # Prepare for run
    ########################
    supply = [supply_ammount] * len(sectors)

    reader = list(opendrift_model.readers.items())[0][1]

    if not isinstance(duration, timedelta):
        duration = timedelta(seconds=duration)

    if time is None:
        time = reader.start_time
    if not isinstance(time, list):
        time = [time]
    if len(time) == 2 and time[1] > time[0] + duration:
        ValueError("Seed time can't be longer than runtime")

    if not time_step_output:
        time_step_output = time_step

    ############################################
    # Seeding element for simulation
    ############################################

    # seed elements within the polygon sectors
    for i, sector in enumerate(sectors):
        z = depth
        opendrift_model.seed_within_polygon(
            lons=sector[0], lats=sector[1], number=supply[i], time=time, origin_marker=i, z=z)

    # get the original lons and lats of the seed
    initial_lons, initial_lats = opendrift_model.get_lonlats()
    initial_lons, initial_lats = copy.copy(initial_lons)[:, 0], copy.copy(initial_lats)[:, 0]

    ##########################
    # Running the simulation
    ##########################

    # run simulation
    opendrift_model.run(duration=duration, steps=steps, end_time=end_time, time_step=time_step,
                        time_step_output=time_step_output)

    print("simulation run completed")

    # returning the origin location of elements, the amount supplied
    # by each sector and the coordinates of the sectors
    lons, lats = opendrift_model.get_lonlats()
    firstlast = np.ma.notmasked_edges(lons, axis=1)
    index_of_last = firstlast[1][1]
    final_lons = lons.T[index_of_last, range(lons.T.shape[1])]
    final_lats = lats.T[index_of_last, range(lons.T.shape[1])]

    return {'lons': initial_lons, 'lats': initial_lats}, {'lons': final_lons, 'lats': final_lats}, supply, sectors


def run_for_connectivity_domain(opendrift_model, delta=None, domain=None,
                                time=None, steps=None, duration=None, end_time=None, time_step=None,
                                depth=0, supply_coords=[], supply_amount=10,
                                _radius=100, time_step_output=None, reader_landmask=reader_global_landmask.Reader()):
    """
    seed and run OpenDrift model before creation of a connectivity map with sectors from a grid of the domain

    :param opendrift_model: opendrift model for the simulation
    :param delta: size of each sector by degrees
    :param domain: a list [xmin, xmax, ymin, ymax] describing the domain for the simulation and connectivity map
    :param time: datenum or list.
        one element - The time at which particles are seeded/released.
        list with two elements - elements are seeded continuously from start/first to end/last time.
        list with more than two elements - the number of elements is equal to len(time) and are seeded as a time series.
    :param steps: integer, maximum number of steps. End of simulation will be self.start_time + steps*self.time_step
    :param duration: timedelta defining the length of the simulation
    :param end_time: datetime object defining the end of the simulation
    :param time_step: interval between particles updates, in seconds or as timedelta. Default: 3600 seconds (1 hour)
    :param depth: depth at which the elements are seeded
    :param supply_coords: specific coordinates for elements seeding
    :param supply_amount: the amount of elements to supply at each sector
    :param _radius: scalar or array, radius in meters around each lon-lat pair, within which particles will be
    randomly seeded.
    :param time_step_output: Time step at which element properties are stored and eventually written to file.
    Timedelta object or seconds. Default: same as time_step, meaning that all steps are stored
    :param reader_landmask: (optional) reader to check if the elements are on land,

    :return: {'lons': initial_lons, 'lats': initial_lats}, {'lons': final_lons, 'lats': final_lats},
    supply, sectors - the origin location of elements, the amount supplied by each sector and the
    coordinates of the sectors
    """
    ########################
    # Prepare for run
    ########################

    print('Calculating connectivity map')

    reader = list(opendrift_model.readers.items())[0][1]
    if isinstance(reader, pyproj.Proj):
        proj = reader
    elif isinstance(reader, str):
        proj = pyproj.Proj(reader)
    else:
        proj = reader.proj

    if not isinstance(duration, timedelta):
        duration = timedelta(seconds=duration)

    if domain is None:
        xmin, xmax, ymin, ymax = reader.xmin, reader.xmax, reader.ymin, reader.ymax
    else:
        xmin, xmax, ymin, ymax = domain

    len_y = ymax - ymin
    len_x = xmax - xmin

    if time is None:
        time = reader.start_time
    if not isinstance(time, list):
        time = [time]
    if len(time) == 2 and time[1] > time[0] + duration:
        ValueError("Seed time can't be longer than runtime")

    if not time_step_output:
        time_step_output = time_step

    # divide the domain to sectors
    xs = np.linspace(xmin, xmax - delta, int(round((xmax - xmin) / delta)))
    ys = np.linspace(ymin, ymax - delta, int(round((ymax - ymin) / delta)))

    X, Y = np.meshgrid(xs, ys)
    # getting the lons and lats of the each sector (the corner of the sector)
    # corrected by the projection
    lons, lats = proj(X, Y, inverse=True)

    # Initializing a list to contain all the sectors and the number of
    # element seeded in each one
    sectors = []
    supply = []

    # Initializing a list to contain all the sectors that are on land,
    # so no elements will be seeded on them
    sectors_on_land = []

    print('Creating sectors')

    # Create the sectors and the supply list
    sec_num = 0
    for i, lon in enumerate(lons[0, :]):
        for j, lat in enumerate(lats[:, 0]):
            # Creating and adding the sector to the list of sectors
            sectors.append(([lon, lon, lon + delta, lon + delta],
                            [lat, lat + delta, lat + delta, lat]))
            supply.append(0)
            # Check if the sector is in sea or on land, and if on land than
            # add it to the sectors_on_land list
            if (reader_landmask.__on_land__(np.array([lon]), np.array([lat])) 
                    and reader_landmask.__on_land__(np.array([lon]), np.array([lat + delta]))
                        and reader_landmask.__on_land__(np.array([lon + delta]), np.array([lat + delta]))
                            and reader_landmask.__on_land__(np.array([lon + delta]), np.array([lat]))):
                sectors_on_land.append(sec_num)
            sec_num += 1

    # If coordinates for supply were given, seed elements only in sectors
    # containing a supply coordinate
    if supply_coords:
        for coords in supply_coords:
            lon = coords[0]
            lat = coords[1]
            # Check if the coordinate is within the domain
            if not (xmax > lon > xmin
                    and ymax > lat > ymin):
                continue
            # find the sector the contains the coordinate
            # this assumes the domain is quadrangular (xmin, xmax, ymin, ymax)
            # and creates a grid where each cell is a sector, than translates the
            # cell index ([grid_x, grid_y]) to a number 's', an index of a one
            # dimensional array
            grid_x_index = np.floor((lon - xmin) / delta)
            grid_y_index = np.floor((lat - ymin) / delta)

            final_sector = round(grid_x_index * (len_y / delta) + grid_y_index)

            try:
                supply[final_sector] += supply_amount
            except:
                print("error:", grid_x_index, grid_y_index, lon,
                      lat, final_sector, "outside of domain")

    # If the seeding should be performed over the entire domain
    # (in all the sectors), add supply to the sectors that are in sea
    else:
        for i in range(len(supply)):
            if i not in sectors_on_land:
                supply[i] += supply_amount

    print('Seeding elements')

    ############################################
    # Seeding element for simulation
    ############################################

    # loop through sectors and seed elements
    for i, sector in enumerate(sectors):
        # check if elements need to be seeded in the sector
        if supply[i]:
            N = supply[i]
            z = depth
            # z = depth*np.random.uniform(0, 1, N)
            opendrift_model.seed_elements(lon=np.average(sector[0]), lat=np.average(sector[1]),
                                          number=N, radius=_radius, time=time, origin_marker=i, z=z)

    # get the origin lons and lats of the seed
    initial_lons, initial_lats = opendrift_model.get_lonlats()
    initial_lons, initial_lats = copy.copy(initial_lons)[:, 0], copy.copy(initial_lats)[:, 0]

    print('Seed complete, running simulation')
    print('Running simulation')

    ##########################
    # Running the simulation
    ##########################

    # run simulation
    opendrift_model.run(duration=duration, steps=steps, end_time=end_time, time_step=time_step,
                        time_step_output=time_step_output)

    print('simulation run completed')

    # returning the origin location of elements, the amount supplied
    # by each sector and the coordinates of the sectors
    lons, lats = opendrift_model.get_lonlats()
    firstlast = np.ma.notmasked_edges(lons, axis=1)
    index_of_last = firstlast[1][1]
    final_lons = lons.T[index_of_last, range(lons.T.shape[1])]
    final_lats = lats.T[index_of_last, range(lons.T.shape[1])]

    return {'lons': initial_lons, 'lats': initial_lats}, {'lons': final_lons, 'lats': final_lats}, supply, sectors


# Connectivity matrix and categories 

def calculate_connectivity_matrix(initial_lonlats, final_lonlats, supply, sectors,
                                  last_z=None, max_depth=None, min_depth=None,
                                  polygons=False):
    """
    Calculating the connectivity map and matrix after running the simulation
    :param first_lonlats: the origin location of the elements (at seed time, provided by run_for_connectivity_domain())
    :param last_lonlats: the location of the elements at the end of the simulation
    :param supply: list of the amount supplied by each sector (provided by run_for_connectivity_domain())
    :param sectors: list of all the sectors, provided by run_for_connectivity_domain())
    :param max_depth: (optional) the max depth of elements for the connectivity map
    :param min_depth: (optional) the min depth of elements for the connectivity map
    :param polygons: is the connectivity from a domain or polygons?

    :return: C, supply - Connectivity matrix, supply for each sector
    """

    #########################
    # Preparation
    #########################

    # check if max_depth or min_depth is provided for depth analysis
    check_depth = False
    if max_depth or min_depth:
        check_depth = True
        if min_depth is None:
            min_depth = max_depth
        elif max_depth is None:
            max_depth = min_depth

    # get the current lons and lats of the elements (after the simulation run)
    initial_lons = initial_lonlats['lons']
    initial_lats = initial_lonlats['lats']
    final_lons = final_lonlats['lons']
    final_lats = final_lonlats['lats']

    # create a copy of the supply array to make the changes locally
    supply_copy = copy.copy(supply)

    print("starting analysis")

    num_sectors = len(sectors)

    # Initializing the connectivity matrix
    # column C[:,i] represents the elements settled in the sector i
    # row C[i,:] represents the elements supplied by the sector i
    #         1_rec   2_rec   ...
    # 1_sup |       |       | ...
    # 2_sup |       |       | ...
    # ...
    C = np.zeros((int(num_sectors), int(num_sectors)))

    ##################################
    # Create connectivity matrix
    ##################################

    # if the sectors are predetermined polygons
    if polygons:
        # create bounding boxes for the polygons for faster intersection check
        boxes = [(np.min(sector[0]), np.max(sector[0]), np.min(sector[1]), np.max(sector[1])) for sector in sectors]
        # loop through the elements to count number of elements in each sector and create connectivity matrix
        for _id in range(len(initial_lons)):
            # get the first and last location of the element
            initial_lon = initial_lons[_id]
            initial_lat = initial_lats[_id]
            final_lon = final_lons[_id]
            final_lat = final_lats[_id]
            if check_depth:
                last_depth = last_z[_id]

            initial_sector, final_sector = -1, -1

            # loop through the sectors and check if the element is inside one of the sectors
            for sector_index, sector in enumerate(sectors):
                # find the origin sector of the element
                if initial_sector == -1 and boxes[sector_index][1] >= initial_lon >= boxes[sector_index][0] and boxes[sector_index][3] >= initial_lat >= boxes[sector_index][2]:
                    if is_point_in_area(sector[0], sector[1], (initial_lon, initial_lat)):
                        initial_sector = sector_index
                # find the sector the element settled in
                if final_sector == -1 and boxes[sector_index][1] >= final_lon >= boxes[sector_index][0] and boxes[sector_index][3] >= final_lat >= boxes[sector_index][2]:
                    if is_point_in_area(sector[0], sector[1], (final_lon, final_lat)):
                        final_sector = sector_index
                # if both sectors were found, stop the check
                if initial_sector != -1 and final_sector != -1:
                    break
            # if the elements settled outside the sectors, ignore the element
            else:
                supply_copy[initial_sector] -= 1
                continue
            # check if the element depth is outside of the selected depth
            if check_depth and not max_depth >= last_depth >= min_depth:
                supply_copy[initial_sector] -= 1
                continue
            # if both sectors were found, add the pair to the connectivity matrix
            if final_sector != -1:
                C[initial_sector, final_sector] += 1

        # return the connectivity matrix
        return C

    delta = round(np.abs(sectors[0][0][1] - sectors[0][0][2]), len(str(sectors[0][0][2]).split('.')[1]))

    # data about the domain
    xmin, xmax, ymin, ymax = sectors[0][0][0], sectors[-1][0][0] + delta, sectors[0][1][0], sectors[-1][1][0] + delta
    n_ymax = ymax - ymin
    n_xmax = xmax - xmin

    # loop through elements to count number of elements in each sector and create connectivity matrix
    for _id in range(len(initial_lons)):
        # get the first and last location of the element
        initial_lon = initial_lons[_id]
        initial_lat = initial_lat[_id]
        final_lon = final_lon[_id]
        final_lat = final_lats[_id]
        if check_depth:
            last_depth = last_z[_id]

        # find the origin sector of the element
        # this assumes the domain is quadrangular (xmin, xmax, ymin, ymax)
        # and creates a grid where each cell is a sector, than translates the
        # cell index ([grid_x, grid_y]) to a number 'i', an index of a one
        # dimensional array
        grid_x = np.floor((initial_lon - xmin) / delta)
        grid_y = np.floor((initial_lat - ymin) / delta)

        # the origin sector
        initial_sector = round(grid_x * (n_ymax / delta) + grid_y)

        # find the sector the element settled in
        grid_x = np.floor((final_lon - xmin) / delta)
        grid_y = np.floor((final_lat - ymin) / delta)

        # the sector the element settled in
        final_sector = round(grid_x * (n_ymax / delta) + grid_y)

        # find elements that were seeded on land or that are out of domain, and remove them from the calculation
        # on land, i.e. the location haven't changed during the simulation
        if final_lon == initial_lon and final_lat == initial_lat:
            supply_copy[initial_sector] -= 1
            continue
        # outside of domain
        if not (xmax > final_lon > xmin
                and ymax > final_lat > ymin):
            supply_copy[initial_sector] -= 1
            continue
        # outside of chosen depth range
        if (check_depth and not
        max_depth >= last_depth >= min_depth):
            supply_copy[initial_sector] -= 1
            continue

        # if both sectors were found, add the pair to the connectivity matrix
        try:
            C[initial_sector, final_sector] += 1
        except:
            print("error:", grid_x, grid_y, final_lon, final_lat, final_sector)
            raise (ValueError)

    print("Connectivity matrix created")

    # return the connectivity matrix
    return C


def categorize_sectors(sectors, C, cat_threshold=5):
    """
    Categorizes the sectors by their settling / supply values.
    This values are based on the paper "Connectivity of larval stages of
    sedentary marine communities between hard substrates and offshore
    structures in the North Sea" by Johan van der Molen, Luz María García-García,
    Paul Whomersley, Alexander Callaway, Paulette E. Posen & Kieran Hyde
    https://www.nature.com/articles/s41598-018-32912-2
    :param sectors: the sectors and their coordinates
    :param cat_threshold: threshold for connectivity categories
    :param C: Connectivity matrix from calculate_connectivity_matrix_from_domain()
    :return: Nsettle_values, Nsupply_values, R_values, Srec_values, Ssup_values, Sector_categories
    """

    num_sectors = len(sectors)

    Nsettle_values = []  # the number of particles settling in each sector
    Nsupply_values = []  # the number of particles supplied by each sector
    R_values = []  # settling / supply factor, Nsupply - Nsettle / Nsupply + Nsettle
    Srec_values = []  # the number of sectors from which each sector received
    Ssup_values = []  # the number of sectors to which each sector supplied
    Sector_categories = []  # connectivity categories (1-8)

    # calculate return values
    for i in range(num_sectors):
        if not i % 100:
            print(i, end=" ")
        # count the number of elements settled in sector i
        n_settle = (C[:, i]).sum()
        Nsettle_values.append(n_settle)
        # count the number of elements supplied by sector i
        n_supply = (C[i, :]).sum()
        Nsupply_values.append(n_supply)
        # calculate the settling / supply factor of sector i
        r = (n_supply - n_settle) / (n_supply + n_settle)
        R_values.append(r)
        # count the number of sectors from which a sector received
        srec = (C[:, i] > 0).sum()  # columns represent the target sector
        Srec_values.append(srec)
        # count the number of sectors to which a sector supplied
        ssup = (C[i, :] > 0).sum()  # rows represent the source sector
        Ssup_values.append(ssup)

        # classify sectors into connectivity categories
        if r > 0.5:
            Sector_categories.append(
                8) if ssup > cat_threshold else Sector_categories.append(7)
        elif -0.5 <= r <= 0.5:
            if ssup > cat_threshold:
                Sector_categories.append(
                    6) if srec > cat_threshold else Sector_categories.append(5)
            else:
                Sector_categories.append(
                    4) if srec > cat_threshold else Sector_categories.append(3)
        elif r < 0.5:
            Sector_categories.append(
                2) if srec > cat_threshold else Sector_categories.append(1)
        else:
            Sector_categories.append(np.nan)

    return Nsettle_values, Nsupply_values, R_values, Srec_values, Ssup_values, Sector_categories


# Plotting 


def plot_connectivity_map(sector_categories=[], sectors=[], corners=[],
                          figsize=(10, 6), title='', nan_color='white', bad_alpha=0,
                          coastline_color='black', ocean_color="white",
                          colormap='viridis',
                          vmin=None, vmax=None, lut=None):
    """
    plots the connectivity categories of each sector

    :param sector_categories:  connectivity categories of the sectors
    :param corners: the corner coordinates for the plot
    :param figsize: tuple, (width, height) in inches.
    :param title: title for the plot
    :param nan_color: color for nan values
    :param bad_alpha: transparency of nan_color
    :param coastline_color: color for the coastline
    :param ocean_color: color for the ocean
    :param land_color: color for the land
    :param colormap: the colormap for the plot
    :param vmin: minimum value of sector_categories
    :param vmax: maximum value of sector_categories
    :param lut: if not None the colormap will be resampled to have lut entries in the lookup table
    """
    # data about the domain
    delta = round(np.abs(sectors[0][0][1] - sectors[0][0][2]), len(str(sectors[0][0][2]).split('.')[1]))
    xmin, xmax, ymin, ymax = sectors[0][0][0], sectors[-1][0][0] + delta, sectors[0][1][0], sectors[-1][1][0] + delta

    if not corners:
        corners = [xmin, xmax, ymin, ymax]

    # divide the map into sectors for plotting
    lon = np.linspace(xmin, xmax - delta, int(round((xmax - xmin) / delta)))
    lat = np.linspace(ymin, ymax - delta, int(round((ymax - ymin) / delta)))

    # set the colormap
    cmap = plt.cm.get_cmap(colormap, lut)
    # set color of nan values
    cmap.set_bad(color=nan_color, alpha=bad_alpha)
    cmap.set_under(color=nan_color, alpha=bad_alpha)

    # projection
    map_proj = cartopy.crs.PlateCarree()

    # plot base map
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=map_proj)

    # reshape the categories for plotting
    sc_ = np.reshape(sector_categories, (-1, len(lat)))

    # set min-max values for color bar
    if not vmin:
        vmin = np.nanmin(sector_categories)
        if vmin == 0:
            vmin = 0.00000001
    if not vmax:
        vmax = np.nanmax(sector_categories)

    # plot the sectors
    topo_plot = ax.pcolormesh((lon + delta / 2),
                              (lat + delta / 2),
                              sc_.T, vmin=vmin, vmax=vmax,
                              transform=map_proj,
                              cmap=cmap)

    ax.coastlines(color=coastline_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=(ocean_color))
    ax.add_feature(cartopy.feature.LAND, edgecolor='black')

    # plot only within the corners
    ax.set_extent(corners)
    # plot the gridlines
    ax.gridlines(draw_labels=True, dms=True,
                 x_inline=False, y_inline=False)

    # add colorbar
    axpos = ax.get_position()
    cbar_ax = fig.add_axes([axpos.x1 + 0.07, axpos.y0, 0.03, axpos.height])
    cbar = fig.colorbar(topo_plot, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    ax.gridlines(draw_labels=True, dms=True,
                 x_inline=False, y_inline=False)

    # set the title
    ax.set_title(title, fontsize=24)
    plt.show()


def plot_connectivity_polygons(sectors=[], sector_categories=[], corners=[],
                               figsize=(10, 6), title='', nan_color='white', bad_alpha=0,
                               coastline_color='black', ocean_color='white', colormap="viridis",
                               vmin=None, vmax=None, lut=None, numbering=False):
    """
    plots the connectivity categories of each sector

    :param sectors: the sectors to plot
    :param sector_categories:  connectivity categories of the sectors
    :param corners: the corner coordinates for the plot
    :param title: title for the plot
    :param nan_color: color for nan values
    :param bad_alpha: transparency of nan_color
    :param coastline_color: color for the coastline
    :param ocean_color: color for the ocean
    :param land_color: color for the land
    :param colormap: the colormap for the plot
    :param vmin: minimum value of sector_categories
    :param vmax: maximum value of sector_categories
    :param lut: if not None the colormap will be resampled to have lut entries in the lookup table
    :param numbering: boolean, if True sector numbers will be annotated over the sectors
    """
    # data about the domain
    if not corners:
        boxes = np.array(
            [[np.min(sector[0]), np.max(sector[0]), np.min(sector[1]), np.max(sector[1])] for sector in sectors])
        corners = [min(boxes[:, 0]), max(boxes[:, 1]), min(boxes[:, 2]), max(boxes[:, 3])]

    # set the colormap
    cmap = plt.get_cmap(colormap, lut)
    # set color of nan values
    cmap.set_bad(color=nan_color, alpha=bad_alpha)

    # projection
    map_proj = cartopy.crs.PlateCarree()

    # plot base map
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=map_proj)
    ax.coastlines(color=coastline_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=(ocean_color))
    ax.add_feature(cartopy.feature.LAND, edgecolor='black')

    # plot only within corners
    ax.set_extent(corners)
    # plot gridlines
    ax.gridlines(draw_labels=True, dms=True,
                 x_inline=False, y_inline=False)

    # plot the sectors from polygons
    for n, sector in enumerate(sectors):
        pgon = Polygon(tuple([(sector[0][i], sector[1][i]) for i in range(len(sector[0]))]))
        ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor=cmap((sector_categories[n]) / 8), alpha=0.8)
        if numbering:
            cx = pgon.representative_point().x
            cy = pgon.representative_point().y
            ax.annotate(str(n), (cx, cy))
        

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 8))
    sm._A = []
    plt.colorbar(sm, ax=ax)



    # add title
    ax.set_title(title, fontsize=24)
    plt.show()

def plot_connectivity_network(sectors=[], connectivity_matrix=[], sector_categories=[], corners=[], 
                               figsize=(10, 6), title='', nan_color='white', bad_alpha=0,
                               coastline_color='black', ocean_color='white', colormap="viridis",
                               vmin=None, vmax=None, lut=None, numbering=False):
    """
    plots the a network between the sectors

    :param sectors: the sectors to plot
    :param connectivity_matrix: the connectivity matrix of the sectors
    :param sector_categories:  connectivity categories of the sectors
    :param corners: the corner coordinates for the plot
    :param title: title for the plot
    :param nan_color: color for nan values
    :param bad_alpha: transparency of nan_color
    :param coastline_color: color for the coastline
    :param ocean_color: color for the ocean
    :param land_color: color for the land
    :param colormap: the colormap for the plot
    :param vmin: minimum value of sector_categories
    :param vmax: maximum value of sector_categories
    :param lut: if not None the colormap will be resampled to have lut entries in the lookup table
    :param numbering: boolean, if True sector numbers will be annotated over the sectors
    """
    # data about the domain
    if not corners:
        boxes = np.array(
            [[np.min(sector[0]), np.max(sector[0]), np.min(sector[1]), np.max(sector[1])] for sector in sectors])
        corners = [min(boxes[:, 0]), max(boxes[:, 1]), min(boxes[:, 2]), max(boxes[:, 3])]

    # set the colormap
    cmap = plt.get_cmap(colormap, lut)
    # set color of nan values
    cmap.set_bad(color=nan_color, alpha=bad_alpha)

    # projection
    map_proj = cartopy.crs.PlateCarree()

    # plot base map
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=map_proj)
    ax.coastlines(color=coastline_color)
    ax.add_feature(cartopy.feature.OCEAN, facecolor=(ocean_color))
    ax.add_feature(cartopy.feature.LAND, edgecolor='black')

    # plot only within corners
    ax.set_extent(corners)
    # plot gridlines
    ax.gridlines(draw_labels=True, dms=True,
                 x_inline=False, y_inline=False)

    # plot the sectors from polygons
    for n, sector in enumerate(sectors):
        pgon = Polygon(tuple([(sector[0][i], sector[1][i]) for i in range(len(sector[0]))]))
        ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor=cmap((sector_categories[n]) / 8), alpha=0.8)
        cx = pgon.representative_point().x
        cy = pgon.representative_point().y
        if numbering:
            ax.annotate(str(n), (cx, cy))
        for tar_n, tar_sector in enumerate(sectors):
            pgon2 = Polygon(tuple([(tar_sector[0][i], tar_sector[1][i]) for i in range(len(tar_sector[0]))]))
            cx2 = pgon2.representative_point().x
            cy2 = pgon2.representative_point().y
            weight = connectivity_matrix[n, tar_n]
            if weight == 0:
                continue
            plt.plot([cx, cx2], [cy, cy2], linewidth=1, color='black', transform= cartopy.crs.Geodetic())

        

    # add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 8))
    sm._A = []
    plt.colorbar(sm, ax=ax)

    # add title
    ax.set_title(title, fontsize=24)
    plt.show()


# for plotting connectivity from polygons using grid
# def plot_connectivity_polygons_(sectors=[], domain=[], sector_categories=[], corners=[],
#                                 figsize=(10, 6), title='', bad_color='white', bad_alpha=0, coastline_color='black'):
#     if corners == []:
#         corners = domain
#     delta = 0.05

#     lon = np.arange(domain[0], domain[1], delta)
#     lat = np.arange(domain[2], domain[3], delta)

#     categs = np.zeros((len(lon), len(lat)))
#     categs.fill(np.nan)

#     for i, _lon in enumerate(lon):
#         for j, _lat in enumerate(lat):
#             for n, sector in enumerate(sectors):
#                 if is_point_in_area(sector[0], sector[1], (_lon, _lat)):
#                     categs[i, j] = sector_categories[n]

#     plot_connectivity_map(domain=domain, delta=delta, sector_categories=categs, corners=corners,
#                           figsize=figsize, title=title, nan_color=bad_color, bad_alpha=bad_alpha,
#                           coastline_color=coastline_color)
