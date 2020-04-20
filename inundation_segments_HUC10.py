#!/usr/bin/env python3

import grass.script as gs
# from grass.pygrass.modules import Module
from grass.pygrass.modules.shortcuts import general as g


def main():
    import os
    # Module("r.lake")
    # Module("r.patch")

    # we want to run this repetetively without deleted the created files
    # os.environ['GRASS_OVERWRITE'] = '1'

    stage_heights = 'stage_heights'
    HAND = 'HAND_clip'
    catchhuc = 'catchhuc'

    # list.files()
    # stage_DT = fread("stage_heights_Sept13-Sept19_2019_12PM.csv", na.strings = "NA")
    # CatchID = stage_DT$CatchID

    # Creating vector of names for "r.lake" output layers for each of the COMID's
    # inundName = paste("i_", CatchID, sep = "")
    g.message("Making inundation layer names...")
    # CatchID = []
    # CatchID = gs.run_command("v.db.select", map=stage_heights, columns="COMID")
    # lenCatchID = len(CatchID)

    CatchID = gs.read_command("v.db.select", flags="c", map=stage_heights, columns="COMID")
    CatchID = (CatchID.split("\n"))
    CatchID = CatchID[0:(len(CatchID)-1)]
    CatchID = [int(i) for i in CatchID]
    print(CatchID)
    lenCatchID = len(CatchID) - 1

    i_string = 'i'
    start_str = list(''.join([i_string] * len(CatchID)))
    # print(start_str)
    inundName = list(map(lambda x, y: x + '_' +  str(y), start_str, CatchID))
    print(inundName)

    g.message("Making v.extract layer names...")
    # Creating vector of names for "v.extract" output layers for each of the COMID's
    # v_COMID = paste("COMID_", CatchID, sep = "")
    COM_string = 'C'
    start_str = list(''.join([COM_string] * len(CatchID)))
    # print(start_str)
    v_COMID = list(map(lambda x, y: x + '_' +  str(y), start_str, CatchID))
    print(v_COMID)

    g.message("Extracting stage heights...")
    # Extracting the stage heights from the shapefile for September 13th
    # stage_13 = []
    # stage_13 = gs.run_command("v.db.select", map=stage_heights, columns="stag_13")
    stage_19 = gs.read_command("v.db.select", flags="c", map=stage_heights, columns="stag_19")
    # stage_13 = (stage_13.split("\r"))
    stage_19 = (stage_19.split("\r\n"))
    stage_19 = stage_19[0:(len(stage_19)-1)]
    print(stage_19)
    lenStage19 = len(stage_19)-1
    # stage_13 = [float(i) for i in stage_13]
    for i in range(0,lenStage19):
        if stage_19[i] == '':
            stage_19[i] = '0.0'


    for i in stage_19:
        if '.' in i:
            float(i)
    print(stage_19)

    COMID_seg = "COMID_seg"
    buff_rast = "buff_rast"

    g.message("Starting for loop...")
    for i in range(0,lenCatchID):
      # gs.read_command("g.region", raster=HAND, nsres=10, ewres=10, n=22070, s=-261550, w=1354140, e=1733390, rows=28362, cols=37925, flags="a")

# north:      22070
# south:      -261550
# west:       1354140
# east:       1733390
# nsres:      10
# ewres:      10
# rows:       28362
# cols:       37925

      tmp_CatchID = CatchID[i]
      # where_exp = paste("COMID == ", tmp_CatchID, sep = "")
      where_exp = "COMID == " + str(tmp_CatchID)

      gs.read_command("v.extract", input=stage_heights, output=v_COMID[i], type="line", where=where_exp, overwrite=True)
      g.message("Finished extraction...starting vector conversion to raster")


      # r.buffer input=COMID out=test_buff_rast dist=3000 units=meters --overwrite
      # r.mask -r
      # r.mask raster=test_buff_rast

      gs.read_command("v.to.rast", input=v_COMID[i], output=COMID_seg, type="line", use="attr", attribute_column="COMID", overwrite=True)
      # gs.read_command("r.mask", flags="r")
      # gs.read_command("r.buffer", input=COMID_seg, out=buff_rast, distances=8000, units="meters")
      # gs.read_command("r.mask", raster=buff_rast)
      # gs.read_command("g.region", raster=buff_rast, nsres=10, ewres=10, flags="a")
      # g.message("Finished vector conversion, buffer,mask, and and chaning g.region...")
      g.message("Finished vector conversion...")

      g.message("Creating mask for subwatershed...")
      gs.run_command('r.mask', raster=catchhuc, maskcats=tmp_CatchID)

      # stage_tmp = stage_DT[CatchID == tmp_CatchID,]
      stage_tmp = stage_19[i]
      g.message("...starting inundation simulation")


      # execGRASS("v.to.rast", input=v_COMID[i], output="stage_13", type="line", use="attr", attribute_column="stage_13", flags=c("overwrite"))
      gs.read_command("r.lake", elevation=HAND, lake=inundName[i], water_level=stage_tmp, seed=COMID_seg, overwrite=True)
      gs.run_command('r.mask', flags='r')
      g.message("Finished inundation simulation, removed mask...now onto the next reach ID")


    g.message("Finished inundation simulation for all segments...patching segments")
    # gs.read_command("r.patch", input=inundName, output="Sept16Flood")
    gs.read_command("r.patch", input=inundName[:-1], output="Sept19Flood", overwrite=True)



    # generate list of flood names: inundNames = inundNames1, inundNames2...
    # for x in length(stage_Oct07$COMID)
        # WL = Stage height @ COMID
        # r.lake(elevation=CyberGIS_HAND water_level=WL lake=inundNames1 seed=flowlines_rast)
    # r.patch(input=inundNames output=Oct_07_12pm_Inund)



if __name__ == '__main__':
    main()
