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
        
        param9 = arcpy.Parameter(
            displayName="Output Raster",
            name="output_raster",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Output"
        )

        return [param0, param1, param2, param3, param4, param5, param6, param7, param8, param9]

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

    def calc_degs(self, dtm, dsm, distance, curvature, vertical_offset):
        return atan(((dsm - curvature) - (dtm + vertical_offset)) / distance) * (180 / pi)

    def get_name(self, path):
        return os.path.splitext(os.path.split(path)[1])[0]

    def deg_to_valid_mils(self, degs):
        m = degs * (6400 / 360)
        if (m > 1600) or (m < -1600):
            return None
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
        dtm = parameters[3].valueAsText

        # Get spatial reference from input surface raster
        sr = arcpy.Describe(parameters[2].valueAsText).spatialReference
        arcpy.env.outputCoordinateSystem = sr

        ao_selection = parameters[0].valueAsText

        if ao_selection == "By View Extent":
            ao_extent = self.get_view_extent()
        else:
            ao_poly = arcpy.Describe(parameters[1].valueAsText).extent
            ao_extent = Extent(ao_poly.XMin, ao_poly.YMin, ao_poly.XMax, ao_poly.YMax)
            
        # Generate fishnet and get center points
        arcpy.AddMessage(f"Using extent: {ao_extent.lowerLeft} | {ao_extent.upperRight}")
        arcpy.SetProgressor("default", "Downsampling deployment area extent...")
        arcpy.CreateFishnet_management(
            r"memory\fn",
            f"{ao_extent.XMin} {ao_extent.YMin}",
            f"{ao_extent.XMin} {ao_extent.YMax}",
            cellsize,
            cellsize,
            corner_coord=f"{ao_extent.XMax} {ao_extent.YMax}"
        )

        ao_cell_centers = r"memory\fn_label"
        arcpy.AddFields_management(
            ao_cell_centers,
            [["distance", "LONG"],
             ["bearing", "LONG"]])
        
        # dtm_sampled will be the final table from which the output raster will be built
        dtm_sampled = Sample(
            dtm,
            ao_cell_centers,
            r"memory\dtm_samp",
            unique_id_field="OID",
            generate_feature_class="FEATURE_CLASS"
        )
        
        # Add the fields required to dtm_sampled
        arcpy.AddFields_management(
            dtm_sampled,
            [["incline_deg", "DOUBLE"],
             ["incline_mil", "DOUBLE"]]
        )
        
        """
        dtm_sampled row now looks like:
        
        (OID, (SHAPE@XY), LOCATIONID, X, Y, <dtm_name>, incline_deg, incline_mil)
        """

        arcpy.SetProgressor("default", "Updating fields...")
        with arcpy.da.UpdateCursor(ao_cell_centers, ["distance", "bearing"]) as cursor:
            for row in cursor:
                row[0] = distance
                row[1] = bearing
                cursor.updateRow(row)
        del cursor

        arcpy.AddXY_management(ao_cell_centers)

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

        """
        Generate points along lines of bearing.  We'll write this out because I've
        found even on very well resourced workstations the script has problems properly
        completing the surface raster sampling from the points in memory.
        """
        arcpy.SetProgressor("default", "Creating sample points along lines...")
        samplepoints = arcpy.GeneratePointsAlongLines_management(
            az_lines,
            os.path.join(scratch, "samplepoints"),
            "DISTANCE",
            Distance=f"{interval} meters",
            Include_End_Points="END_POINTS"
        )
        arcpy.AddXY_management(samplepoints)

        # Sample the DSM
        arcpy.SetProgressor("default", "Sampling surface dataset...")
        dsm_sampled = Sample(dsm, samplepoints, os.path.join(scratch, "dsm_sam"),
                             generate_feature_class="FEATURE_CLASS")

        # Add from_origin, dsm, and curvature to samplepoints
        arcpy.SetProgressor("default", "Adding analysis fields to sampled points...")
        arcpy.AddFields_management(
            samplepoints,
            [["from_origin", "FLOAT"],
             ["dsm", "FLOAT"],
             ["curvature", "DOUBLE"]]
        )

        """
        Update each point's distance from origin as a member of its ORIG_FID group.  This will
        be used to calculate curvature value.

        Copy over the DSM sample values to samplepoints.
        """
        arcpy.SetProgressor("default", "Calculating each point's curvature and distance from its origin...")
        c = defaultdict(list)
        with arcpy.da.SearchCursor(samplepoints, ["OID@", "ORIG_FID", "from_origin"]) as sc:
            for row in sc:
                c[row[1]].append(row)
        c = json.loads(json.dumps(c))
        del sc

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
        del uc

        # Copy over DSM sampled values to samplepoints
        self.copy_sample_values(dsm_sampled, samplepoints, "dsm")

        """
        Schema of samplepoints now looks like:

        (OID, (SHAPE@XY), ORIG_FID, POINT_X, POINT_Y, DISTANCE, BEARING, SHAPE_LENGTH, FROM_ORIGIN, DSM, CURVATURE)
        """
        # Calculate the inclination of each point in samplepoints relative to its origin in dtm_sampled
        
        s = defaultdict(list)
        with arcpy.da.SearchCursor(samplepoints, ["ORIG_FID", "from_origin", "dsm", "curvature"]) as sp_cursor:
            for row in sp_cursor:
                s[row[0]].append(row)
        del sp_cursor
        
        samps = dict((k, [i for i in v if i[1] > 0]) for k, v in s.items())
        
        with arcpy.da.UpdateCursor(dtm_sampled, "*") as dtm_cursor:
            for row in dtm_cursor:
                loc_id = row[2]
                dtm_val = row[5]
                if loc_id in samps.keys():
                    i = []
                    s_points = samps[loc_id]
                    for sp in s_points:
                        d = sp[1]
                        c = sp[3]
                        dsm_val = sp[2]
                        try:
                            inclination = self.calc_degs(dtm_val, dsm_val, d, c, vert_offset)
                            i.append(inclination)
                        except TypeError:  # Use no-data value for distances not covered by surface raster
                            i.append(-9999)
                    max_inclination = max(i)
                    row[6] = max_inclination
                    row[7] = self.deg_to_valid_mils(max_inclination)
                dtm_cursor.updateRow(row)
        del dtm_cursor

        # Now create the final output raster from the points based on inclination
        arcpy.SetProgressor("default", "Writing inclination raster...")
        if parameters[9].valueAsText:
            out_raster = parameters[9].valueAsText
        else:
            out_raster = os.path.join(default_db, arcpy.CreateScratchName(
                prefix=f"Incl_{bearing}_",
                suffix="",
                data_type="RasterDataset"
            ))
        
        incl_raster = arcpy.PointToRaster_conversion(
            dtm_sampled, 
            "incline_mil", 
            out_raster, 
            cell_assignment="MEAN",
            cellsize=cellsize
        )

        return incl_raster
