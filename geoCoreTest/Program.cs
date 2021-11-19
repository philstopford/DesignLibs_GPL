using geoCoreLib;
using geoLib;
using System;
using System.Collections.Generic;

namespace geoCoreTest;

internal class Program
{
    private static string baseDir = "D:\\Google Drive\\Semi\\geocore_test\\";

    private static void Main(string[] args)
    {
        //test_basic();
        //test_cellrefarray_basic();
        // test_cellrefarray_nested();
        test_cell_export();
        test_cell_export_complex();
    }

    private static void test_cell_export()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray\\";

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();

        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
    }

    private static void test_cell_export_complex()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray\\";

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested_2.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();

        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        List<GCPolygon> t3a = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(layer: 1, datatype: 0);
        List<GCPolygon> t3b = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(layer: 2, datatype: 0);

        List<GCPolygon> t4 = gcGDS.convertToPolygons();
        List<GCPolygon> t4a = gcGDS.convertToPolygons(activeLDOnly: true);
        List<GCPolygon> t4b = gcGDS.convertToPolygons(layer: 2, datatype: 0);

        if (t3.Count != t4.Count)
        {
            throw new Exception("test_cell_export_complex failed : t3,t4");
        }

        if (t3a.Count != t4a.Count)
        {
            throw new Exception("test_cell_export_complex failed : t3a, t4a");
        }

        if (t3b.Count != t4b.Count)
        {
            throw new Exception("test_cell_export_complex failed : t3b, t4b");
        }
    }

    private static void test_cellrefarray_basic()
    {
        string arrayDir = baseDir + "cellrefarray\\";

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // The array is in cell 'a'
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GeoLibPointF[]> geo = gcGDS.points(flatten: true);
    }

    private static void test_cellrefarray_nested()
    {
        string arrayDir = baseDir + "cellrefarray\\";

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested_2.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // The array is in cell 'a'
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GeoLibPointF[]> geo2 = gcGDS.points(flatten: true);
    }

    private static void test_1()
    {

        string gdsDir = baseDir + "gdsFiles\\";
        string oasDir = baseDir + "oasFiles\\";
        string outDir = baseDir + "geoCoreTestOut\\";

        // Basic tests.
        defineAndWrite_Box(outDir);
        defineAndWrite_Polygon(outDir);
        defineAndWrite_Cellref(outDir);
        defineAndWrite_CellrefArray(outDir);
        defineAndWrite_Circle(outDir);
        defineAndWrite_Path(outDir);
        defineAndWrite_Text(outDir);
        defineAndWrite_CTrapezoid(outDir);
        defineAndWrite_Trapezoid(outDir);

        defineAndWrite_Box_LayerDataType(1, 0, outDir);
        defineAndWrite_Box_LayerDataType(128, 0, outDir);
        defineAndWrite_Box_LayerDataType(255, 0, outDir);
        defineAndWrite_Box_LayerDataType(256, 0, outDir);
        defineAndWrite_Box_LayerDataType(511, 0, outDir);
        defineAndWrite_Box_LayerDataType(512, 0, outDir);
        defineAndWrite_Box_LayerDataType(16383, 0, outDir);
        defineAndWrite_Box_LayerDataType(16384, 0, outDir);
        defineAndWrite_Box_LayerDataType(32767, 0, outDir);
        defineAndWrite_Box_LayerDataType(32768, 0, outDir);
        defineAndWrite_Box_LayerDataType(1, 128, outDir);
        defineAndWrite_Box_LayerDataType(1, 255, outDir);
        defineAndWrite_Box_LayerDataType(1, 256, outDir);
        defineAndWrite_Box_LayerDataType(1, 511, outDir);
        defineAndWrite_Box_LayerDataType(1, 512, outDir);
        defineAndWrite_Box_LayerDataType(1, 32767, outDir);
        defineAndWrite_Box_LayerDataType(1, 32768, outDir);

        defineAndWrite_manyCellrefs(4, outDir);

        string baseName = "paisley";
        loadSaveTest(baseName, gdsDir, oasDir, outDir);
    }

    private static void loadSaveTest(string baseName, string gdsDir, string oasDir, string outDir)
    {
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsDir + baseName + ".gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // Can we write consistent with what was read.
        gds.gdsWriter gw = new(gcGDS, outDir + baseName + "_out.gds");
        gw.save();

        // Check that we can write out the GDS-sourced drawing to Oasis, proving the internal drawing is handled consistently.
        oasis.oasWriter gow = new(gcGDS, outDir + baseName + "_GDSout.oas");
        gow.save();

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasDir + baseName + ".oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();

        // Can we write consistent with what was read.
        oasis.oasWriter ow = new(gcOAS, outDir + baseName + "_out.oas");
        ow.save();

        // Check that we can write out the OAS-sourced drawing to GDS, proving the internal drawing is handled consistently.
        gds.gdsWriter ogw = new(gcOAS, outDir + baseName + "_OASout.gds");
        ogw.save();
    }

    private static void defineAndWrite_Box(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 1E-3, // 0.001 / 1E-6;
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addBox(0, 0, 10, 20, 1, 0);

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_box.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_box.oas");
        ow.save();
    }

    private static void defineAndWrite_Polygon(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        // L
        GeoLibPoint[] poly = new GeoLibPoint[6];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(0, 20);
        poly[2] = new GeoLibPoint(10, 20);
        poly[3] = new GeoLibPoint(10, 10);
        poly[4] = new GeoLibPoint(20, 10);
        poly[5] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // triangle
        poly = new GeoLibPoint[3];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(10, 20);
        poly[2] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 2, 0);

        // pentagram
        poly = new GeoLibPoint[5];
        poly[0] = new GeoLibPoint(5, 0);
        poly[1] = new GeoLibPoint(0, 10);
        poly[2] = new GeoLibPoint(10, 20);
        poly[3] = new GeoLibPoint(20, 10);
        poly[4] = new GeoLibPoint(15, 0);

        gcell.addPolygon(poly, 3, 0);

        // trapezoid
        poly = new GeoLibPoint[4];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(5, 20);
        poly[2] = new GeoLibPoint(15, 20);
        poly[3] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 4, 0);

        // parallelogram
        poly = new GeoLibPoint[4];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(10, 20);
        poly[2] = new GeoLibPoint(20, 20);
        poly[3] = new GeoLibPoint(10, 0);

        gcell.addPolygon(poly, 5, 0);

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_polygon.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_polygon.oas");
        ow.save();
    }

    private static void defineAndWrite_Path(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        GeoLibPoint[] path = new GeoLibPoint[5];
        path[0] = new GeoLibPoint(0, 0);
        path[1] = new GeoLibPoint(0, 10);
        path[2] = new GeoLibPoint(20, 10);
        path[3] = new GeoLibPoint(20, 40);
        path[4] = new GeoLibPoint(0, 40);

        gcell.addPath(path, 1, 0);
        gcell.elementList[gcell.elementList.Count - 1].setWidth(5); // note that Oasis only supports factors of 2, so this gets rounded down to make a 4 unit path at the writer (for Oasis).
        gcell.elementList[gcell.elementList.Count - 1].setCap(2); // extend the path by half the width at the line ends.

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_path.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_path.oas");
        ow.save();
    }

    private static void defineAndWrite_Circle(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addCircle(1, 0, new GeoLibPoint(10, 10), 5.0);

        g.setDrawing(drawing_);
        g.setValid(true);

        // GDS looks janky for circles; no special category unlike the Oasis format. Snaps to grid and looks ugly.
        gds.gdsWriter gw = new(g, outDir + "simple_circle.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_circle.oas");
        ow.save();
    }

    private static void defineAndWrite_Cellref(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        GeoLibPoint[] poly = new GeoLibPoint[6];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(0, 20);
        poly[2] = new GeoLibPoint(10, 20);
        poly[3] = new GeoLibPoint(10, 10);
        poly[4] = new GeoLibPoint(20, 10);
        poly[5] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 1, 0);

        bool mirror_x = false;
        gcell = drawing_.addCell();
        gcell.addCellref();

        gcell.elementList[gcell.elementList.Count - 1].setPos(new GeoLibPoint(10, 0));
        gcell.elementList[gcell.elementList.Count - 1].setCellRef(drawing_.findCell("test"));
        gcell.elementList[gcell.elementList.Count - 1].setName("test");
        gcell.elementList[gcell.elementList.Count - 1].rotate(0);
        gcell.elementList[gcell.elementList.Count - 1].scale(2);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[gcell.elementList.Count - 1].setMirrorx();
                break;
        }

        gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test2";

        poly = new GeoLibPoint[6];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(0, 30);
        poly[2] = new GeoLibPoint(10, 30);
        poly[3] = new GeoLibPoint(10, 10);
        poly[4] = new GeoLibPoint(20, 10);
        poly[5] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 1, 0);

        mirror_x = true;
        gcell = drawing_.addCell();
        gcell.addCellref();
        gcell.elementList[gcell.elementList.Count - 1].setPos(new GeoLibPoint(20, 20));
        gcell.elementList[gcell.elementList.Count - 1].setCellRef(drawing_.findCell("test2"));
        gcell.elementList[gcell.elementList.Count - 1].setName("test2");
        gcell.elementList[gcell.elementList.Count - 1].rotate(0);
        gcell.elementList[gcell.elementList.Count - 1].scale(2);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[gcell.elementList.Count - 1].setMirrorx();
                break;
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_cellref.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_cellref.oas");
        ow.save();
    }

    private static void defineAndWrite_CellrefArray(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        GeoLibPoint[] poly = new GeoLibPoint[6];
        poly[0] = new GeoLibPoint(0, 0);
        poly[1] = new GeoLibPoint(0, 20);
        poly[2] = new GeoLibPoint(10, 20);
        poly[3] = new GeoLibPoint(10, 10);
        poly[4] = new GeoLibPoint(20, 10);
        poly[5] = new GeoLibPoint(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        GeoLibPoint[] array = new GeoLibPoint[3];
        array[0] = new GeoLibPoint(0, 0);
        array[1] = new GeoLibPoint(100, 0);
        array[2] = new GeoLibPoint(0, 80);

        bool mirror_x = false;
        gcell = drawing_.addCell();
        gcell.addCellref(drawing_.findCell("test"), new GeoLibPoint(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[gcell.elementList.Count - 1].setPos(new GeoLibPoint(0, 0));
        gcell.elementList[gcell.elementList.Count - 1].setName("test");
        gcell.elementList[gcell.elementList.Count - 1].rotate(0);
        gcell.elementList[gcell.elementList.Count - 1].scale(1);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[gcell.elementList.Count - 1].setMirrorx();
                break;
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_cellrefarray.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_cellrefarray.oas");
        ow.save();
    }

    private static void defineAndWrite_Text(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addText(1, 0, new GeoLibPoint(10, 10), "Text");

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_text.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_text.oas");
        ow.save();
    }

    private static void defineAndWrite_CTrapezoid(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        int w = 40;
        int h = 20;
        int x = 10;
        int y = 20;

        for (int i = 0; i < 25; i++)
        {
            GeoLibPoint[] pa = new GeoLibPoint[5];

            int[,] coords = i switch
            {
                0 => new int[4, 4]
                {
                    {0, 0, 0, 0}, // x=0*w+0*h, y=0*w+0*h ...
                    {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}
                },
                1 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                2 => new int[4, 4] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                3 => new int[4, 4] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                4 => new int[4, 4] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
                5 => new int[4, 4] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                6 => new int[4, 4] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                7 => new int[4, 4] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
                8 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
                9 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                10 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
                11 => new int[4, 4] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                12 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 1, 0}},
                13 => new int[4, 4] {{0, 0, 1, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                14 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
                15 => new int[4, 4] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
                16 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
                17 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
                18 => new int[4, 4] {{0, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
                19 => new int[4, 4] {{0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}},
                20 => new int[4, 4] {{0, 0, 0, 0}, {0, 1, 0, 1}, {0, 2, 0, 0}, {0, 0, 0, 0}},
                21 => new int[4, 4] {{0, 0, 0, 1}, {0, 2, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 1}},
                22 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 2, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
                23 => new int[4, 4] {{1, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 2, 0}, {1, 0, 0, 0}},
                24 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                25 => new int[4, 4] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}},
                _ => new int[4, 4]
            };

            for (int pt = 0; pt < 4; pt++)
            {
                x = 0;
                if (coords[pt, 0] != 0)
                {
                    x += coords[pt, 0] * w;
                }

                if (coords[pt, 1] != 0)
                {
                    x += coords[pt, 1] * h;
                }

                y = 0;
                if (coords[pt, 2] != 0)
                {
                    y += coords[pt, 2] * w;
                }

                if (coords[pt, 3] != 0)
                {
                    y += coords[pt, 3] * h;
                }

                pa[pt] = new GeoLibPoint(x, y);

                if (x > w)
                {
                    w = x;
                }
                if (y > h)
                {
                    h = y;
                }
            }

            pa[pa.Length - 1] = new GeoLibPoint(pa[0]);

            gcell.addPolygon(pa, i + 1, 0);
        }


        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_ctrapezoid.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_ctrapezoid.oas");
        ow.save();
    }

    private static void defineAndWrite_Trapezoid(string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test"; int x = 10;
        int y = 20;
        int w = 40;
        int h = 20;

        bool trapezoid_orientation = false;
        int trapezoid_delta_a = 5;
        int trapezoid_delta_b = 5;

        GeoLibPoint[] pa = new GeoLibPoint[5];

        switch (trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new GeoLibPoint(x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new GeoLibPoint(x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new GeoLibPoint(x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new GeoLibPoint(x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new GeoLibPoint(x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new GeoLibPoint(x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new GeoLibPoint(x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new GeoLibPoint(x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new GeoLibPoint(pa[0]);

        gcell.addPolygon(pa, 1, 0);

        trapezoid_orientation = true;

        switch (trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new GeoLibPoint(x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new GeoLibPoint(x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new GeoLibPoint(x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new GeoLibPoint(x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new GeoLibPoint(x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new GeoLibPoint(x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new GeoLibPoint(x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new GeoLibPoint(x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new GeoLibPoint(pa[0]);

        gcell.addPolygon(pa, 2, 0);

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "simple_trapezoid.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "simple_trapezoid.oas");
        ow.save();
    }

    private static void defineAndWrite_Box_LayerDataType(int layer, int datatype, string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 1E-3, // 0.001 / 1E-6;
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addBox(0, 0, 10, 20, layer, datatype);

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + "L" + layer + "D" + datatype + "_box.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + "L" + layer + "D" + datatype + "_box.oas");
        ow.save();
    }

    private static void defineAndWrite_manyCellrefs(int edge, string outDir)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 0.001,
            libname = "noname"
        };

        bool mirror_x = false;
        GCCell master_gcell = drawing_.addCell();

        for (int i = 0; i < edge * edge; i++)
        {
            GCCell gcell = drawing_.addCell();
            gcell.accyear = 2018;
            gcell.accmonth = 12;
            gcell.accday = 5;
            gcell.acchour = 2;
            gcell.accmin = 10;
            gcell.accsec = 10;
            gcell.modyear = 2018;
            gcell.modmonth = 12;
            gcell.modday = 5;
            gcell.modhour = 2;
            gcell.modmin = 10;
            gcell.modsec = 10;

            gcell.cellName = "test" + i;

            GeoLibPoint[] poly = new GeoLibPoint[6];
            poly[0] = new GeoLibPoint(0, 0);
            poly[1] = new GeoLibPoint(0, 20);
            poly[2] = new GeoLibPoint(10, 20);
            poly[3] = new GeoLibPoint(10, 10);
            poly[4] = new GeoLibPoint(20, 10);
            poly[5] = new GeoLibPoint(20, 0);

            gcell.addPolygon(poly, 1, 0);

            master_gcell.addCellref();
            GCElement cellref = master_gcell.elementList[master_gcell.elementList.Count - 1];
            cellref.setPos(new GeoLibPoint(40 * (i % edge), 40 * Math.Floor((double)i / edge)));
            cellref.setCellRef(drawing_.findCell("test" + i));
            cellref.setName("test" + i);
            cellref.rotate(0);
            cellref.scale(1);
            switch (mirror_x)
            {
                case true:
                    cellref.setMirrorx();
                    break;
            }
        }
        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, outDir + edge + "_cellref.gds");
        gw.save();

        oasis.oasWriter ow = new(g, outDir + edge + "_cellref.oas");
        ow.save();
    }
}