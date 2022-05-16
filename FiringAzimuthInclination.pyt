import arcpy
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import CreateConstantRaster, Extent, Sample
from collections import defaultdict
import json
from math import *
import os


class Toolbox(object):
    def __init__(self):
        self.label = "DIVARTY Tools"
        self.alias = ""
        self.tools = [FiringAzimuthInclination]


class FiringAzimuthInclination(object):

    def __init__(self):
        """
        Derives angle of inclination over a distance/azimuth based on surface and terrain
        models
        """
        self.name = "FiringAzimuthInclination"
        self.label = "Calcluate Firing Azimuth Inclination"
        self.alias = ""
        self.description = "Calculate inclination by azimuth and distance"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """
        Define parameters
        :return: list
        """
        param0 = arcpy.Parameter(
            displayName="Area of Operations",
            name="area_of_operations",
            datatype="GPString",
            parameterType="Required",
            direction="Input"
        )
        param0.filter.type = "ValueList"
        param0.filter.list = ["By Polygon", "By View Extent"]
        param0.value = "By View Extent"

        # Optional parameter to select an extent based on a polygon
        param1 = arcpy.Parameter(
            displayName="AO Layer",
            name="ao_layer",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input"
        )
        param1.filter.list = ["Polygon"]

        param2 = arcpy.Parameter(
            displayName="Surface Raster",
            name="surface_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input"
        )

        param3 = arcpy.Parameter(
            displayName="Terrain Raster",
            name="terrain_raster",
            datatype="GPRasterLayer",
            parameterType="Required",
            direction="Input"
        )

        param4 = arcpy.Parameter(
            displayName="Cell Size",
            name="cell_size",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )

        param5 = arcpy.Parameter(
            displayName="Distance",
            name="distance",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )

        param6 = arcpy.Parameter(
            displayName="Bearing",
            name="bearing",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )

        param7 = arcpy.Parameter(
            displayName="Interval",
            name="interval",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )
        param7.value = 20  # Sane default

        param8 = arcpy.Parameter(
            displayName="Vertical Offset",
            name="vertical_offset",
            datatype="GPLong",
            parameterType="Required",
            direction="Input"
        )
        param8.value = 2  # Sane default

        return [param0, param1, param2, param3, param4, param5, param6, param7, param8]

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        if parameters[0].valueAsText == "By View Extent":
            parameters[1].enabled = False
            for p in parameters[2:]:
                p.enabled = True
        else:
            for p in parameters[1:]:
                p.enabled = True
        return True

    def updateMessages(self, parameters):
        return True

    def get_view_extent(self):
        p = arcpy.mp.ArcGISProject("CURRENT")
        ext = p.activeView.camera.getExtent()
        ll = ext.lowerLeft
        ur = ext.upperRight
        return Extent(ll.X, ll.Y, ur.X, ur.Y)

    def copy_sample_values(self, source_lyr, target_lyr, target_field):
        """
        Helper method to avoid costly joins.  Extract join key and sample value to dict,
        then call da.UpdateCursor() to update target point feature class with sampled
        raster values.
        """
        sample_vals = dict([(r[2], r[-1]) for r in arcpy.da.SearchCursor(source_lyr, "*")])
        with arcpy.da.UpdateCursor(target_lyr, ["OID@", target_field]) as update_rows:
            for row in update_rows:
                row[1] = sample_vals[row[0]]
                update_rows.updateRow(row)
        del update_rows
        return True

    def curvature(self, sample_distance: int):
        return 0.2032 * ((sample_distance / 1609) ** 2)

    def calc_degs(self, distance, curvature, dsm_val, dtm_val, vertical_offset):
        rads = atan(((dsm_val - curvature) - (dtm_val + vertical_offset)) / distance) * (180 / pi)
        return degrees(rads)

    def get_name(self, path):
        return os.path.splitext(os.path.split(path)[1])[0]

    def deg_to_valid_mils(self, degs):
        m = degs * (6400 / 360)
        if m > 1600:
            return 1600
        elif m < -1600:
            return -1600
        else:
            return m

    def execute(self, parameters, messages):
        # Environments
        arcpy.env.overwriteOutput = True
        p = arcpy.mp.ArcGISProject("CURRENT")
        default_db = p.defaultGeodatabase
        scratch = arcpy.env.scratchGDB

        # Parameters
        cellsize = parameters[4].value
        distance = parameters[5].value
        bearing = parameters[6].value
        interval = parameters[7].value
        vert_offset = parameters[8].value

        dsm = parameters[2].valueAsText
        dsm_name = self.get_name(dsm)
        dtm = parameters[3].valueAsText
        dtm_name = self.get_name(dtm)

        # Get spatial reference from input surface raster
        sr = arcpy.Describe(parameters[2].valueAsText).spatialReference
        arcpy.env.outputCoordinateSystem = sr

        ao_selection = parameters[0].valueAsText

        if ao_selection == "By View Extent":
            ao_extent = self.get_view_extent()
        else:
            ao_poly = arcpy.Describe(parameters[1].valueAsText).extent
            ao_extent = Extent(ao_poly.XMin, ao_poly.YMin, ao_poly.XMax, ao_poly.YMax)

        # Create constant raster denoting area of operations (ao) from ao_extent
        arcpy.SetProgressor("default", "Building constant raster...")
        ao = CreateConstantRaster(1, "FLOAT", cell_size=cellsize, extent=ao_extent)
        ao.save(os.path.join(scratch, "cr"))

        """
        We are doing five things in this block:

        1. Get cell centers as feature class,
        2. Add distance and bearin fields,
        3. Populate each row's distance and bearing,
        4. Add XY coordinates to attributes,
        5. Convert it to a table for deriving the lines of bearing
        """
        arcpy.SetProgressor("default", "Converting raster to points...")
        ao_cell_centers = arcpy.RasterToPoint_conversion(ao, os.path.join(scratch, "rcc"))
        arcpy.AddFields_management(
            ao_cell_centers,
            [["distance", "LONG"],
             ["bearing", "LONG"]])

        arcpy.SetProgressor("default", "Adding fields...")
        with arcpy.da.UpdateCursor(ao_cell_centers, ["distance", "bearing"]) as cursor:
            for row in cursor:
                row[0] = distance
                row[1] = bearing
                cursor.updateRow(row)

        arcpy.AddXY_management(ao_cell_centers)

        arcpy.SetProgressor("default", "Converting points to table...")
        celltable = arcpy.TableToTable_conversion(ao_cell_centers, r"memory", "celltable")

        # Now construct lines of bearing from the table
        arcpy.SetProgressor("default", "Building lines of bearing...")
        az_lines = arcpy.BearingDistanceToLine_management(
            celltable,
            r"memory\az_lines",
            "POINT_X",
            "POINT_Y",
            distance_field="distance",
            distance_units="METERS",
            bearing_field="bearing",
            spatial_reference=sr
        )

        # Build points along line at fixed interval
        arcpy.SetProgressor("default", "Creating sample points along lines...")
        samplepoints = arcpy.GeneratePointsAlongLines_management(
            az_lines,
            r"memory\samplepoints",
            "DISTANCE",
            Distance=f"{interval} meters",
            Include_End_Points="END_POINTS"
        )

        # Sample the DSM
        arcpy.SetProgressor("default", "Sampling surface dataset...")
        dsm_sampled = Sample(dsm, samplepoints, r"memory\dsm_sam", generate_feature_class="FEATURE_CLASS")

        # Sample the DTM
        arcpy.SetProgressor("default", "Sampling bare earth dataset...")
        dtm_sampled = Sample(dtm, samplepoints, r"memory\dtm_sam", generate_feature_class="FEATURE_CLASS")

        # Add from_origin, dsm, dtm, curvature, incl_deg and incl_mils to samplepoints
        arcpy.SetProgressor("default", "Adding analysis fields to sampled points...")
        arcpy.AddFields_management(
            samplepoints,
            [["from_origin", "FLOAT"],
             ["dsm", "FLOAT"],
             ["dtm", "FLOAT"],
             ["curvature", "DOUBLE"],
             ["incl_deg", "DOUBLE"],
             ["incl_mils", "DOUBLE"]])

        """
        Update each point's distance from origin as a member of its ORIG_FID group.  This will
        be used to calculate curvature value.

        Copy over the DSM/DTM sample values to samplepoints.
        """
        arcpy.SetProgressor("default", "Calculating each point's curvature and distance from its origin...")
        c = defaultdict(list)
        with arcpy.da.SearchCursor(samplepoints, ["OID@", "ORIG_FID", "from_origin"]) as sc:
            for row in sc:
                c[row[1]].append(row)
        c = json.loads(json.dumps(c))

        # Calculate distance to origin for each point
        sp_with_intervals = {int(k): [(a, b, interval * i) for i, (a, b, c) in enumerate(v)] for k, v in c.items()}

        with arcpy.da.UpdateCursor(samplepoints, ["OID@", "ORIG_FID", "from_origin", "curvature"]) as uc:
            for u_row in uc:
                if u_row[1] in sp_with_intervals.keys():
                    for n_row in sp_with_intervals[u_row[1]]:
                        if n_row[0] == u_row[0]:  # Check OID alignment
                            u_row[2] = n_row[2]
                            u_row[3] = self.curvature(sample_distance=u_row[2])
                uc.updateRow(u_row)

        # Copy over DSM/DTM sampled values to samplepoints
        self.copy_sample_values(dsm_sampled, samplepoints, "dsm")
        self.copy_sample_values(dtm_sampled, samplepoints, "dtm")

        """
        Schema of samplepoints now looks like:

        (OID, (SHAPE@XY), ORIG_FID, POINT_X, POINT_Y, DISTANCE, BEARING, SHAPE_LENGTH, FROM_ORIGIN, DSM, DTM, CURVATURE, INCL_DEG, INCL_MILS)
        """
        # Calculate inclination of each sampled point
        arcpy.SetProgressor("default", "Calculating inclination in degrees and mils...")
        with arcpy.da.UpdateCursor(samplepoints, "*") as uc:
            for row in uc:
                d = row[8]  # Distance from origin point
                c = row[11]  # Curvature
                dsm_val = row[9]  # Sampled DSM value
                dtm_val = row[10]  # Sampled DTM value
                try:
                    if d:
                        if d == 0:
                            row[12] = float(0)  # This is the origin point, set incl_deg to 0
                            row[13] = float(0)  # This is the origin point, set incl_mils to 0
                            uc.updateRow(row)
                        else:
                            row[12] = self.calc_degs(d, c, dsm_val, dtm_val, vert_offset)
                            row[13] = self.deg_to_valid_mils(row[12])
                            uc.updateRow(row)
                    else:
                        # d should always be populated, otherwise something went very wrong upstream
                        raise TypeError
                except Exception as te:
                    # Remove try/except for production, replace with more graceful error handling
                    arcpy.AddWarning(te)
                    pass

        # Now create the final output raster from the points based on inclination
        arcpy.SetProgressor("default", "Writing inclination raster...")
        out_raster_path = os.path.join(default_db,
                                       arcpy.CreateScratchName(prefix="FAzIncl_", suffix="", data_type="RasterDataset"))
        incl_raster = arcpy.PointToRaster_conversion(samplepoints, "incl_mils", out_raster_path, cell_assignment="MEAN",
                                                     cellsize=cellsize)

        return incl_raster
