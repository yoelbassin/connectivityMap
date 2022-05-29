from __future__ import annotations
from typing import List
import typing
import numpy as np
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon

Point = typing.NamedTuple('Point', [('lon', float), ('lat', float)])

class Sector:
    def __init__(self, lons: List[float], lats: List[float]) -> None:
        if not len(lons) == len(lats):
            raise ValueError(
                "expected lons and lats to be the same size, got {} lons and {} lats".format(
                    len(lons), len(lats)
                )
            )
        self._lons = lons
        self._lats = lats
        self._vertex_num = len(self._lons)

        self._n_settled = None
        self._n_supplied = None
        self._r = None
        self._s_received = None
        self._s_supplied = None
        self._category = None

    @classmethod
    def create_sectors_domain(cls, domain: List[int], delta: float) -> List[Sector]:
        """
        creates a list of Sector objects from a rectangular domain, with squares of size delta

        :param domain: A list containing the 4 corners of the domain
        :param delta: The size of each rectangle in the domain

        :return: List[Sector] - list of rectangular sectors
        """
        xmin, xmax, ymin, ymax = domain

        # divide the domain to sectors
        xs = np.linspace(xmin, xmax - delta, int(round(abs(xmax - xmin)) / delta))
        ys = np.linspace(ymin, ymax - delta, int(round(abs(ymax - ymin)) / delta))

        X, Y = np.meshgrid(xs, ys)
        # # getting the lons and lats of the each sector (the corner of the sector)
        # # corrected by the projection
        # lons, lats = proj(X, Y, inverse=True)
        lons, lats = X, Y

        # Initializing a list to contain all the sectors and the number of
        # element seeded in each one
        sectors: List[Sector] = []

        print('Creating sectors')

        # Create the sectors and the supply list
        for _, lon in enumerate(lons[0, :]):
            for _, lat in enumerate(lats[:, 0]):
                # Creating and adding the sector to the list of sectors
                sectors.append(cls([lon, lon, lon + delta, lon + delta],
                                [lat, lat + delta, lat + delta, lat]))

        return sectors

    @classmethod
    def sectors_from_shapefile(cls, shapefile: str) -> List[Sector]:
        """
        creates a list of Sector objects from a shapefile

        :param shapefile: str - Shapefile path

        :return: List[Sector] list of polygon sectors
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

                r = geom.GetGeometryRef(0)
                lons = [r.GetY(j) for j in range(r.GetPointCount())]
                lats = [r.GetX(j) for j in range(r.GetPointCount())]
                if not lons or not lats:
                    continue
                sectors.append(cls(lons, lats))
        return sectors

    def check_if_inside(self, point: Point) -> bool:
        """
        Checks if a point is inside the sector

        :param point: (float, float)

        :return: bool
        """
        if type(point) != tuple or len(point) != 2:
            raise ValueError("A vertex is a tuple of length 2")
        nvert = len(self._lons)
        vertx = self._lons
        verty = self._lats
        testx, testy = point[0], point[1]
        c = 0
        j = nvert - 1
        for i in range(0, nvert):
            if ((verty[i] > testy) != (verty[j] > testy)) and (
                testx
                < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i])
                + vertx[i]
            ):
                c = not (c)
            j = i
        return c

    def set_values(
        self,
        n_settled: int,
        n_suplied: int,
        r: float,
        s_received: int,
        s_supplied: int,
        category: float,
    ) -> None:
        self._n_settled = n_settled
        self._n_supplied = n_suplied
        self._r = r
        self._s_received = s_received
        self._s_supplied = s_supplied
        self._category = category

    @property
    def lons(self) -> List[float]:
        return self._lons

    @property
    def lats(self) -> List[float]:
        return self._lats

    @property
    def n_settled(self) -> int:
        return self._n_settled

    @property
    def n_supplied(self) -> int:
        return self._n_supplied

    @property
    def r(self) -> float:
        return self._r

    @property
    def s_received(self) -> int:
        return self._s_received

    @property
    def s_supplied(self) -> int:
        return self._s_supplied

    @property
    def category(self) -> int:
        return self._category

    def __getitem__(self, num: int) -> List[float]:
        if num == 0:
            return self.lons
        if num == 1:
            return self.lats
        raise ValueError("expected a 1 or 2, got {}".format({num}))


class Connectivity:
    def __init__(
        self, lons: List[float], lats: List[float], sectors: List[Sector], cat_threshold: int = 5
    ) -> None:
        """
        Create a connectivity object

        For trajectories - We only need the starting and finising points.

        :param lons: List[float] - A list containing items longitude trajectories
        :param lats: List[float] - A list containing items latitude trajectories
        :param cat_threshold: int (default 5) - Threshold for connectivity categories
        """
        firstlast = np.ma.notmasked_edges(lons, axis=1)
        index_of_last = firstlast[1][1]
        final_lons = lons.T[index_of_last, range(lons.T.shape[1])]
        final_lats = lats.T[index_of_last, range(lons.T.shape[1])]

        self._lons = list(zip(lons[:, 0], final_lons))
        self._lats = list(zip(lats[:, 0], final_lats))
        self._sectors = sectors

        self._calc()
        self._categorize_sectors(cat_threshold)

    def _calc(self):
        """
        Create connectivity matrix for the sectors
        """

        num_sectors = len(self._sectors)

        # Initializing the connectivity matrix
        # column C[:,i] represents the elements settled in the sector i
        # row C[i,:] represents the elements supplied by the sector i
        #         1_rec   2_rec   ...
        # 1_sup |       |       | ...
        # 2_sup |       |       | ...
        # ...
        connectivity_matrix = np.zeros((int(num_sectors), int(num_sectors)))

        # create bounding boxes for the polygons for faster intersection check
        boxes = [
            (
                np.min(sector.lons),
                np.max(sector.lons),
                np.min(sector.lats),
                np.max(sector.lats),
            )
            for sector in self._sectors
        ]
        # loop through the elements to count number of elements in each sector and create connectivity matrix
        for _id in range(len(self._lons)):
            # get the first and last location of the element
            initial_lon, final_lon = self._lons[_id]
            initial_lat, final_lat = self._lats[_id]

            initial_sector, final_sector = -1, -1

            # loop through the sectors and check if the element is inside one of the sectors
            for sector_index, sector in enumerate(self._sectors):
                # find the origin sector of the element
                if (
                    initial_sector == -1
                    and boxes[sector_index][1] >= initial_lon >= boxes[sector_index][0]
                    and boxes[sector_index][3] >= initial_lat >= boxes[sector_index][2]
                ):
                    if sector.check_if_inside((initial_lon, initial_lat)):
                        initial_sector = sector_index
                # find the sector the element settled in
                if (
                    final_sector == -1
                    and boxes[sector_index][1] >= final_lon >= boxes[sector_index][0]
                    and boxes[sector_index][3] >= final_lat >= boxes[sector_index][2]
                ):
                    if sector.check_if_inside((final_lon, final_lat)):
                        final_sector = sector_index
                # if both sectors were found, stop the check
                if initial_sector != -1 and final_sector != -1:
                    break
            # # check if the element depth is outside of the selected depth
            # if check_depth and not max_depth >= last_depth >= min_depth:
            #     supply_copy[initial_sector] -= 1
            #     continue

            # if both sectors were found, add the pair to the connectivity matrix
            if final_sector != -1 and initial_sector != -1:
                connectivity_matrix[initial_sector, final_sector] += 1

        self._connectivity_matrix = connectivity_matrix

    def _categorize_sectors(self, cat_threshold=5):
        """
        Categorizes the sectors by their settling / supply values.
        This values are based on the paper "Connectivity of larval stages of
        sedentary marine communities between hard substrates and offshore
        structures in the North Sea" by Johan van der Molen, Luz María García-García,
        Paul Whomersley, Alexander Callaway, Paulette E. Posen & Kieran Hyde
        https://www.nature.com/articles/s41598-018-32912-2
        
        :param cat_threshold: threshold for connectivity categories
        """

        C = self._connectivity_matrix

        # calculate return values
        for i, sector in enumerate(self._sectors):
            if not i % 100:
                print(i, end=" ")
            # count the number of elements settled in sector i
            n_settle = (C[:, i]).sum()
            # count the number of elements supplied by sector i
            n_supply = (C[i, :]).sum()

            if n_supply + n_settle == 0:
                r = np.nan
            else:
                # calculate the settling / supply factor of sector i
                r = (n_supply - n_settle) / (n_supply + n_settle)
            # count the number of sectors from which a sector received
            srec = (C[:, i] > 0).sum()  # columns represent the target sector
            # count the number of sectors to which a sector supplied
            ssup = (C[i, :] > 0).sum()  # rows represent the source sector

            # classify sectors into connectivity categories
            if r > 0.5:
                if ssup > cat_threshold:
                    category = 8
                else:
                    category = 7
            elif -0.5 <= r <= 0.5:
                if ssup > cat_threshold:
                    if srec > cat_threshold:
                        category = 6
                    else:
                        category = 5
                else:
                    if srec > cat_threshold:
                        category = 4
                    else:
                        category = 3
            elif r < 0.5:
                if srec > cat_threshold:
                    category = 2
                else:
                    category = 1
            else:
                category = np.nan

            sector.set_values(n_settle, n_supply, r, srec, ssup, category)

    def plot(
        self,
        sector_categories: List[float],
        buffer: float =0.1,
        corners: List[float]=None,
        squares: bool=False,
        delta: float=None,
        title: str="",
        nan_color: str="white",
        bad_alpha: float=0,
        coastline_color: str="black",
        land_color: str=cfeature.COLORS["land"],
        ocean_color: str="white",
        colormap: str="viridis",
        vmin: float=None,
        vmax: float=None,
        lut=None,
    ):
        """
        plots the connectivity categories of each sector

        :param sector_categories: List[float] - the value of each sector
        :param buffer: float (dofault 0.1) - margin to add to the corners
        :param corners: List[float] - the corner coordinates for the plot
        :param sqares: bool (default False) - 
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
        if corners:
            lonmin, lonmax, latmin, latmax = corners
        else:
            lonmin = np.nanmin(self._lons)
            lonmax = np.nanmax(self._lons)
            latmin = np.nanmin(self._lats)
            latmax = np.nanmax(self._lats)
            lonmin = lonmin - buffer * 2
            lonmax = lonmax + buffer * 2
            latmin = latmin - buffer
            latmax = latmax + buffer

        crs = ccrs.Mercator()
        lscale = "auto"

        meanlat = (latmin + latmax) / 2
        aspect_ratio = float(latmax - latmin) / (float(lonmax - lonmin))
        aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if aspect_ratio > 1:
            fig = plt.figure(figsize=(11.0 / aspect_ratio, 11.0))
        else:
            fig = plt.figure(figsize=(11.0, 11.0 * aspect_ratio))

        ax = fig.add_subplot(111, projection=crs)  # need '111' for Python 2
        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

        # set the colormap
        cmap = plt.cm.get_cmap(colormap, lut).copy()
        # set color of nan values
        cmap.set_bad(color=nan_color, alpha=bad_alpha)
        cmap.set_under(color=nan_color, alpha=bad_alpha)


        # reshape the categories for plotting

        # set min-max values for color bar
        if not vmin:
            vmin = np.nanmin(sector_categories)
            if vmin == 0:
                vmin = 0.00000001
        if not vmax:
            vmax = np.nanmax(sector_categories)


        if squares:
            if not delta:
                raise ValueError(
                    "if in sqares mode, delta (square size) should be provided"
                )
            # divide the map into sectors for plotting
            lon = np.linspace(
                lonmin, lonmax - delta, int(round(abs(lonmax - lonmin) / delta))
            )
            lat = np.linspace(
                latmin, latmax - delta, int(round(abs(latmax - latmin) / delta))
            )

            sc_ = np.reshape(sector_categories, (-1, len(lat)))

            topo_plot = ax.pcolormesh(
                (lon + delta / 2),
                (lat + delta / 2),
                sc_.T,
                vmin=vmin,
                vmax=vmax,
                transform=cartopy.crs.PlateCarree(),
                cmap=cmap,
            )

            # # add colorbar
            axpos = ax.get_position()
            cbar_ax = fig.add_axes([axpos.x1 + 0.07, axpos.y0, 0.03, axpos.height])
            fig.colorbar(topo_plot, cax=cbar_ax)
        
        else:
            for n, sector in enumerate(self._sectors):
                if sector_categories[n] in [np.nan, 0]:
                    continue
                pgon = Polygon(tuple([(sector[0][i], sector[1][i]) for i in range(len(sector[0]))]))
                ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor=cmap((sector_categories[n]) / max(sector_categories)))
                

            # add colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, max(sector_categories)))
            sm._A = []
            plt.colorbar(sm, ax=ax)

        f = cfeature.GSHHSFeature(scale=lscale, levels=[1], facecolor=land_color)
        ax.add_geometries(
            f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
            ccrs.PlateCarree(),
            facecolor=land_color,
            edgecolor=coastline_color,
        )

        gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True)
        if cartopy.__version__ < "0.18.0":
            gl.xlabels_top = False  # Cartopy < 0.18
        else:
            gl.top_labels = None  # Cartopy >= 0.18

        fig.canvas.draw()
        # fig.set_tight_layout(True)

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        # set the title
        ax.set_title(title, fontsize=24)

        ax.plot(transform=cartopy.crs.Geodetic())

    @property
    def n_settle_values(self) -> List[int]:
        # the number of particles settling in each sector
        return [sector.n_settled for sector in self._sectors]

    @property
    def n_supply_values(self) -> List[int]:
        # the number of particles supplied by each sector
        return [sector.n_supplied for sector in self._sectors]

    @property
    def r_values(self) -> List[int]:
        # settling / supply factor, Nsupply - Nsettle / Nsupply + Nsettle
        return [sector.r for sector in self._sectors]

    @property
    def s_received_values(self) -> List[int]:
        # the number of sectors from which each sector received
        return [sector.s_received for sector in self._sectors]

    @property
    def s_supply_values(self) -> List[int]:
        # the number of sectors to which each sector supplied
        return [sector.s_supplied for sector in self._sectors]

    @property
    def categories(self) -> List[int]:
        # connectivity categories (1-8)
        return [sector.category for sector in self._sectors]