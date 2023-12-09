using System.Buffers.Binary;
using Clipper2Lib;
using geoCoreLib;
using utility;

namespace UnitTests;

public class GeoCoreTests
{
    static string baseDir = "/d/development/DesignLibs_GPL/geocore_test/";
    static string outDir = baseDir + "out/";

    [SetUp]
    public static void GeoCoreSetUp()
    {
        valid_test();
        box_test();
        test_circle();
        test_box();
        consistency_from_oasis();
        consistency_from_gds();

        defineAndWrite_Polygon();
        defineAndWrite_Cellref();
        defineAndWrite_CellrefArray();
        defineAndWrite_Circle();
        defineAndWrite_Path();
        defineAndWrite_Text();
        defineAndWrite_CTrapezoid();
        defineAndWrite_Trapezoid();
        defineAndWrite_manyCellrefs();

        defineAndWrite_Box_LayerDataType();
        
        test_cellrefarray_basic();
        test_cellrefarray_nested();
        test_cell_export();
        test_cell_export_complex();
        
    }

    [Test]
    public static void valid_test()
    {
        GeoCore g = new();
        g.reset();
        if (File.Exists(outDir + "should_fail.gds"))
        {
            File.Delete(outDir + "should_fail.gds");
        }
        if (File.Exists(outDir + "should_fail.oas"))
        {
            File.Delete(outDir + "should_fail.oas");
        }
        // These should fail because the GeoCore set-up is invalid.
        try
        {
            gds.gdsWriter gw = new(g, outDir + "/should_fail.gds");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.AreEqual("Provided GeoCore instance is not marked as valid", e.Message);
        }
        try
        {
            oasis.oasWriter ow = new(g, outDir + "/should_fail.oas");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.AreEqual("Provided GeoCore instance is not marked as valid", e.Message);
        }
        
        GCDrawingfield drawing_ = new("test")
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
        
        g.setDrawing(drawing_);
        // These should fail because the GeoCore set-up is invalid.
        try
        {
            gds.gdsWriter gw = new(g, outDir + "/should_fail.gds");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.AreEqual("Provided GeoCore instance is not marked as valid", e.Message);
        }
        try
        {
            oasis.oasWriter ow = new(g, outDir + "/should_fail.oas");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.AreEqual("Provided GeoCore instance is not marked as valid", e.Message);
        }
    }
    
    [Test]
    public static void box_test()
    {
        int dimension = 10;
        int scale_from_nm = 1; // makes a um scale feature.
        for (int x = -10; x < 20; x += 10)
        {
            for (int y = -10; y < 20; y += 10)
            {
                string filename = "box_" + x + "_" + y + "_" + dimension + "_" + dimension;
                int scale = 100; // for 0.01 nm resolution.
                GeoCore g = new();
                g.reset();
                GCDrawingfield drawing_ = new("test")
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
                    databaseunits = 1000 * scale,
                    userunits = 0.001 / scale,
                    libname = "noname"
                };
                Assert.AreEqual(drawing_.accyear, 2018);

                GCCell gcell = drawing_.addCell();
                // Force the below for comparison sake
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

                gcell.addBox(x * scale_from_nm * scale, y * scale_from_nm * scale, dimension * scale_from_nm * scale, dimension * scale_from_nm * scale, 1, 1);

                g.setDrawing(drawing_);
                g.setValid(true);

                if (File.Exists(outDir + "/" + filename + ".gds"))
                {
                    File.Delete(outDir + "/" + filename + ".gds");
                }

                gds.gdsWriter gw = new(g, outDir + "/" + filename + ".gds");
                gw.save();
                Assert.True(File.Exists(outDir + "/" + filename + ".gds"));

                if (File.Exists(outDir + "/" + filename + ".oas"))
                {
                    File.Delete(outDir + "/" + filename + ".oas");
                }

                oasis.oasWriter ow = new(g, outDir + "/" + filename + ".oas");
                ow.save();
                Assert.True(File.Exists(outDir + "/" + filename + ".oas"));

                // Load the files in to see what we have.

                GeoCoreHandler gH_GDS = new();
                gH_GDS.updateGeoCoreHandler(outDir + "/" + filename + ".gds", GeoCore.fileType.gds);
                GeoCore gcGDS = gH_GDS.getGeo();
                Assert.True(gcGDS.isValid());

                GCDrawingfield drawing_gds = gcGDS.getDrawing();
                drawing_gds.databaseunits = 1000 * scale;
                drawing_gds.userunits = 0.001 / scale;
                GCCell cell_gds = drawing_gds.findCell("test");
                int elementCount = cell_gds.elementList.Count;

                string out_filename = outDir + "/" + filename + "_resave_from_gds.gds";
                save_gdsii(gcGDS, out_filename);

                out_filename = outDir + "/" + filename + "_resave_from_gds.oas";
                save_oasis(gcGDS, out_filename);

                GeoCoreHandler gH_OAS = new();
                gH_OAS.updateGeoCoreHandler(outDir + "/" + filename + ".oas", GeoCore.fileType.oasis);
                GeoCore gcOAS = gH_OAS.getGeo();
                Assert.True(gcOAS.isValid());

                GCDrawingfield drawing_oas = gcOAS.getDrawing();
                drawing_oas.databaseunits = 1000 * scale;
                drawing_oas.userunits = 0.001 / scale;
                GCCell cell_oas = drawing_oas.findCell("test");
                elementCount = cell_oas.elementList.Count;

                out_filename = outDir + "/" + filename + "_resave_from_oas.gds";
                save_gdsii(gcOAS, out_filename);

                out_filename = outDir + "/" + filename + "_resave_from_oas.oas";
                save_oasis(gcOAS, out_filename);
            }
        }
    }

    private static void save_gdsii(GeoCore gc, string out_filename)
    {
        if (File.Exists(out_filename))
        {
            File.Delete(out_filename);
        }
        gds.gdsWriter ow = new(gc, out_filename);
        ow.save();
        Assert.True(File.Exists(out_filename));
    }

    private static void save_oasis(GeoCore gc, string out_filename)
    {
        if (File.Exists(out_filename))
        {
            File.Delete(out_filename);
        }
        oasis.oasWriter ow = new(gc, out_filename);
        ow.save();
        Assert.True(File.Exists(out_filename));
    }
    
    [Test]
    public static void test_circle()
    {
        int numberOfPoints = 60;

        double deg_per_pt = 360.0 / numberOfPoints;

        double radius = 5;

        PathD circleD = new();

        for (int pt = 0; pt < numberOfPoints; pt++)
        {
            double deg = deg_per_pt * pt;
            double rad = Utils.toRadians(deg);

            double x = Math.Cos(rad) * radius;
            double y = Math.Sin(rad) * radius;
            
            circleD.Add(new(x, y));
        }

        int scale = 100; // for 0.01 nm resolution.
        
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
            databaseunits = 1000 * scale,
            userunits = 0.001 / scale,
            libname = "noname"
        };
        Assert.AreEqual(drawing_.accyear, 2018);

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
        
        gcell.addPolygon(Clipper.ScalePath64(circleD, scale*scale), 1, 0);
        // For comparison
        gcell.addCircle(1, 1, new(0,0), 500);

        // Do we have the cell in the drawing?
        Assert.AreEqual(drawing_.findCell("test"), gcell);
       
        g.setDrawing(drawing_);
        g.setValid(true);
        
        if (File.Exists(outDir + "poly10_circle11.gds"))
        {
            File.Delete(outDir + "poly10_circle11.gds");
        }
        gds.gdsWriter gw = new(g, outDir + "poly10_circle11.gds");
        gw.save();
        Assert.True(File.Exists(outDir + "poly10_circle11.gds"));

        if (File.Exists(outDir + "poly10_circle11.oas"))
        {
            File.Delete(outDir + "poly10_circle11.oas");
        }
        oasis.oasWriter ow = new(g, outDir + "poly10_circle11.oas");
        ow.save();
        Assert.True(File.Exists(outDir + "poly10_circle11.oas"));
    }

    [Test]
    public static void test_box()
    {

        int scale = 100; // for 0.01 nm resolution.
        
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
            databaseunits = 1000 * scale,
            userunits = 0.001 / scale,
            libname = "noname"
        };
        Assert.AreEqual(drawing_.accyear, 2018);

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

        PathD boxPoly = Clipper.MakePath(new double[]
        {
            0.081, 0.46,
            0.081, 0.47,
            0.1, 0.47,
            0.1, 0.46,
        });
        
        gcell.addPolygon(Clipper.ScalePath64(boxPoly, 1000*scale), 1, 2);
        
        gcell.addBox(9050, 46500, 1900, 1000, 1, 1);

        // Do we have the cell in the drawing?
        Assert.AreEqual(drawing_.findCell("test"), gcell);
       
        g.setDrawing(drawing_);
        g.setValid(true);
        
        if (File.Exists(outDir + "box.gds"))
        {
            File.Delete(outDir + "box.gds");
        }
        gds.gdsWriter gw = new(g, outDir + "box.gds");
        gw.save();
        Assert.True(File.Exists(outDir + "box.gds"));

        if (File.Exists(outDir + "box.oas"))
        {
            File.Delete(outDir + "box.oas");
        }
        oasis.oasWriter ow = new(g, outDir + "box.oas");
        ow.save();
        Assert.True(File.Exists(outDir + "box.oas"));
    }

    [Test]
    public static void test_cell_export()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");
        Assert.AreEqual(gcGDS.activeStructure, 0);

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");
        Assert.AreEqual(gcGDS.activeLD, 0);

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.AreEqual(t1.Count, 1);
        Assert.AreEqual(t1[0].layer_nr, 1);
        Assert.AreEqual(t1[0].datatype_nr, 0);
        Assert.AreEqual(t1[0].pointarray.Count, 7);
        Assert.AreEqual(t1[0].pointarray[0], new Point64(0,0,0));
        Assert.AreEqual(t1[0].pointarray[1], new Point64(0,100,0));
        Assert.AreEqual(t1[0].pointarray[2], new Point64(40,100,0));
        Assert.AreEqual(t1[0].pointarray[3], new Point64(40,40,0));
        Assert.AreEqual(t1[0].pointarray[4], new Point64(100,40,0));
        Assert.AreEqual(t1[0].pointarray[5], new Point64(100,0,0));
        Assert.AreEqual(t1[0].pointarray[6], new Point64(0,0,0));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.AreEqual(gcGDS.activeStructure, 1);
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.AreEqual(t2.Count, 4);
        Assert.AreEqual(t2[0].layer_nr, 1);
        Assert.AreEqual(t2[0].datatype_nr, 0);
        Assert.AreEqual(t2[0].pointarray.Count, 7);
        Assert.AreEqual(t2[0].pointarray[0], new Point64(0,0,0));
        Assert.AreEqual(t2[0].pointarray[1], new Point64(0,100,0));
        Assert.AreEqual(t2[0].pointarray[2], new Point64(40,100,0));
        Assert.AreEqual(t2[0].pointarray[3], new Point64(40,40,0));
        Assert.AreEqual(t2[0].pointarray[4], new Point64(100,40,0));
        Assert.AreEqual(t2[0].pointarray[5], new Point64(100,0,0));
        Assert.AreEqual(t2[0].pointarray[6], new Point64(0,0,0));
        Assert.AreEqual(t2[1].pointarray.Count, 7);
        Assert.AreEqual(t2[1].pointarray[0], new Point64(0,110,0));
        Assert.AreEqual(t2[1].pointarray[1], new Point64(0,210,0));
        Assert.AreEqual(t2[1].pointarray[2], new Point64(40,210,0));
        Assert.AreEqual(t2[1].pointarray[3], new Point64(40,150,0));
        Assert.AreEqual(t2[1].pointarray[4], new Point64(100,150,0));
        Assert.AreEqual(t2[1].pointarray[5], new Point64(100,110,0));
        Assert.AreEqual(t2[1].pointarray[6], new Point64(0,110,0));
        Assert.AreEqual(t2[2].pointarray.Count, 7);
        Assert.AreEqual(t2[2].pointarray[0], new Point64(110,0,0));
        Assert.AreEqual(t2[2].pointarray[1], new Point64(110,100,0));
        Assert.AreEqual(t2[2].pointarray[2], new Point64(150,100,0));
        Assert.AreEqual(t2[2].pointarray[3], new Point64(150,40,0));
        Assert.AreEqual(t2[2].pointarray[4], new Point64(210,40,0));
        Assert.AreEqual(t2[2].pointarray[5], new Point64(210,0,0));
        Assert.AreEqual(t2[2].pointarray[6], new Point64(110,0,0));
        Assert.AreEqual(t2[3].pointarray.Count, 7);
        Assert.AreEqual(t2[3].pointarray[0], new Point64(110,110,0));
        Assert.AreEqual(t2[3].pointarray[1], new Point64(110,210,0));
        Assert.AreEqual(t2[3].pointarray[2], new Point64(150,210,0));
        Assert.AreEqual(t2[3].pointarray[3], new Point64(150,150,0));
        Assert.AreEqual(t2[3].pointarray[4], new Point64(210,150,0));
        Assert.AreEqual(t2[3].pointarray[5], new Point64(210,110,0));
        Assert.AreEqual(t2[3].pointarray[6], new Point64(110,110,0));
        
        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        Assert.AreEqual(gcGDS.activeStructure, 2);
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.AreEqual(t3.Count, 16);
        Assert.AreEqual(t3[0].layer_nr, 1);
        Assert.AreEqual(t3[0].datatype_nr, 0);
        Assert.AreEqual(t3[0].pointarray.Count, 7);
        Assert.AreEqual(t3[0].pointarray[0], new Point64(0,0,0));
        Assert.AreEqual(t3[0].pointarray[1], new Point64(0,100,0));
        Assert.AreEqual(t3[0].pointarray[2], new Point64(40,100,0));
        Assert.AreEqual(t3[0].pointarray[3], new Point64(40,40,0));
        Assert.AreEqual(t3[0].pointarray[4], new Point64(100,40,0));
        Assert.AreEqual(t3[0].pointarray[5], new Point64(100,0,0));
        Assert.AreEqual(t3[0].pointarray[6], new Point64(0,0,0));
        Assert.AreEqual(t3[1].pointarray.Count, 7);
        Assert.AreEqual(t3[1].pointarray[0], new Point64(0,110,0));
        Assert.AreEqual(t3[1].pointarray[1], new Point64(0,210,0));
        Assert.AreEqual(t3[1].pointarray[2], new Point64(40,210,0));
        Assert.AreEqual(t3[1].pointarray[3], new Point64(40,150,0));
        Assert.AreEqual(t3[1].pointarray[4], new Point64(100,150,0));
        Assert.AreEqual(t3[1].pointarray[5], new Point64(100,110,0));
        Assert.AreEqual(t3[1].pointarray[6], new Point64(0,110,0));
        Assert.AreEqual(t3[2].pointarray.Count, 7);
        Assert.AreEqual(t3[2].pointarray[0], new Point64(110,0,0));
        Assert.AreEqual(t3[2].pointarray[1], new Point64(110,100,0));
        Assert.AreEqual(t3[2].pointarray[2], new Point64(150,100,0));
        Assert.AreEqual(t3[2].pointarray[3], new Point64(150,40,0));
        Assert.AreEqual(t3[2].pointarray[4], new Point64(210,40,0));
        Assert.AreEqual(t3[2].pointarray[5], new Point64(210,0,0));
        Assert.AreEqual(t3[2].pointarray[6], new Point64(110,0,0));
        Assert.AreEqual(t3[3].pointarray.Count, 7);
        Assert.AreEqual(t3[3].pointarray[0], new Point64(110,110,0));
        Assert.AreEqual(t3[3].pointarray[1], new Point64(110,210,0));
        Assert.AreEqual(t3[3].pointarray[2], new Point64(150,210,0));
        Assert.AreEqual(t3[3].pointarray[3], new Point64(150,150,0));
        Assert.AreEqual(t3[3].pointarray[4], new Point64(210,150,0));
        Assert.AreEqual(t3[3].pointarray[5], new Point64(210,110,0));
        Assert.AreEqual(t3[3].pointarray[6], new Point64(110,110,0));
        Assert.AreEqual(t3[4].pointarray.Count, 7);
        Assert.AreEqual(t3[4].pointarray[0], new Point64(0,220,0));
        Assert.AreEqual(t3[4].pointarray[1], new Point64(0,320,0));
        Assert.AreEqual(t3[4].pointarray[2], new Point64(40,320,0));
        Assert.AreEqual(t3[4].pointarray[3], new Point64(40,260,0));
        Assert.AreEqual(t3[4].pointarray[4], new Point64(100,260,0));
        Assert.AreEqual(t3[4].pointarray[5], new Point64(100,220,0));
        Assert.AreEqual(t3[4].pointarray[6], new Point64(0,220,0));
        Assert.AreEqual(t3[5].pointarray.Count, 7);
        Assert.AreEqual(t3[5].pointarray[0], new Point64(0,330,0));
        Assert.AreEqual(t3[5].pointarray[1], new Point64(0,430,0));
        Assert.AreEqual(t3[5].pointarray[2], new Point64(40,430,0));
        Assert.AreEqual(t3[5].pointarray[3], new Point64(40,370,0));
        Assert.AreEqual(t3[5].pointarray[4], new Point64(100,370,0));
        Assert.AreEqual(t3[5].pointarray[5], new Point64(100,330,0));
        Assert.AreEqual(t3[5].pointarray[6], new Point64(0,330,0));
        Assert.AreEqual(t3[6].pointarray.Count, 7);
        Assert.AreEqual(t3[6].pointarray[0], new Point64(110,220,0));
        Assert.AreEqual(t3[6].pointarray[1], new Point64(110,320,0));
        Assert.AreEqual(t3[6].pointarray[2], new Point64(150,320,0));
        Assert.AreEqual(t3[6].pointarray[3], new Point64(150,260,0));
        Assert.AreEqual(t3[6].pointarray[4], new Point64(210,260,0));
        Assert.AreEqual(t3[6].pointarray[5], new Point64(210,220,0));
        Assert.AreEqual(t3[6].pointarray[6], new Point64(110,220,0));
        Assert.AreEqual(t3[7].pointarray.Count, 7);
        Assert.AreEqual(t3[7].pointarray[0], new Point64(110,330,0));
        Assert.AreEqual(t3[7].pointarray[1], new Point64(110,430,0));
        Assert.AreEqual(t3[7].pointarray[2], new Point64(150,430,0));
        Assert.AreEqual(t3[7].pointarray[3], new Point64(150,370,0));
        Assert.AreEqual(t3[7].pointarray[4], new Point64(210,370,0));
        Assert.AreEqual(t3[7].pointarray[5], new Point64(210,330,0));
        Assert.AreEqual(t3[7].pointarray[6], new Point64(110,330,0));
        Assert.AreEqual(t3[8].pointarray.Count, 7);
        Assert.AreEqual(t3[8].pointarray[0], new Point64(220,0,0));
        Assert.AreEqual(t3[8].pointarray[1], new Point64(220,100,0));
        Assert.AreEqual(t3[8].pointarray[2], new Point64(260,100,0));
        Assert.AreEqual(t3[8].pointarray[3], new Point64(260,40,0));
        Assert.AreEqual(t3[8].pointarray[4], new Point64(320,40,0));
        Assert.AreEqual(t3[8].pointarray[5], new Point64(320,0,0));
        Assert.AreEqual(t3[8].pointarray[6], new Point64(220,0,0));
        Assert.AreEqual(t3[9].pointarray.Count, 7);
        Assert.AreEqual(t3[9].pointarray[0], new Point64(220,110,0));
        Assert.AreEqual(t3[9].pointarray[1], new Point64(220,210,0));
        Assert.AreEqual(t3[9].pointarray[2], new Point64(260,210,0));
        Assert.AreEqual(t3[9].pointarray[3], new Point64(260,150,0));
        Assert.AreEqual(t3[9].pointarray[4], new Point64(320,150,0));
        Assert.AreEqual(t3[9].pointarray[5], new Point64(320,110,0));
        Assert.AreEqual(t3[9].pointarray[6], new Point64(220,110,0));
        Assert.AreEqual(t3[10].pointarray.Count, 7);
        Assert.AreEqual(t3[10].pointarray[0], new Point64(330,0,0));
        Assert.AreEqual(t3[10].pointarray[1], new Point64(330,100,0));
        Assert.AreEqual(t3[10].pointarray[2], new Point64(370,100,0));
        Assert.AreEqual(t3[10].pointarray[3], new Point64(370,40,0));
        Assert.AreEqual(t3[10].pointarray[4], new Point64(430,40,0));
        Assert.AreEqual(t3[10].pointarray[5], new Point64(430,0,0));
        Assert.AreEqual(t3[10].pointarray[6], new Point64(330,0,0));
        Assert.AreEqual(t3[11].pointarray.Count, 7);
        Assert.AreEqual(t3[11].pointarray[0], new Point64(330,110,0));
        Assert.AreEqual(t3[11].pointarray[1], new Point64(330,210,0));
        Assert.AreEqual(t3[11].pointarray[2], new Point64(370,210,0));
        Assert.AreEqual(t3[11].pointarray[3], new Point64(370,150,0));
        Assert.AreEqual(t3[11].pointarray[4], new Point64(430,150,0));
        Assert.AreEqual(t3[11].pointarray[5], new Point64(430,110,0));
        Assert.AreEqual(t3[11].pointarray[6], new Point64(330,110,0));
        Assert.AreEqual(t3[12].pointarray.Count, 7);
        Assert.AreEqual(t3[12].pointarray[0], new Point64(220,220,0));
        Assert.AreEqual(t3[12].pointarray[1], new Point64(220,320,0));
        Assert.AreEqual(t3[12].pointarray[2], new Point64(260,320,0));
        Assert.AreEqual(t3[12].pointarray[3], new Point64(260,260,0));
        Assert.AreEqual(t3[12].pointarray[4], new Point64(320,260,0));
        Assert.AreEqual(t3[12].pointarray[5], new Point64(320,220,0));
        Assert.AreEqual(t3[12].pointarray[6], new Point64(220,220,0));
        Assert.AreEqual(t3[13].pointarray.Count, 7);
        Assert.AreEqual(t3[13].pointarray[0], new Point64(220,330,0));
        Assert.AreEqual(t3[13].pointarray[1], new Point64(220,430,0));
        Assert.AreEqual(t3[13].pointarray[2], new Point64(260,430,0));
        Assert.AreEqual(t3[13].pointarray[3], new Point64(260,370,0));
        Assert.AreEqual(t3[13].pointarray[4], new Point64(320,370,0));
        Assert.AreEqual(t3[13].pointarray[5], new Point64(320,330,0));
        Assert.AreEqual(t3[13].pointarray[6], new Point64(220,330,0));
        Assert.AreEqual(t3[14].pointarray.Count, 7);
        Assert.AreEqual(t3[14].pointarray[0], new Point64(330,220,0));
        Assert.AreEqual(t3[14].pointarray[1], new Point64(330,320,0));
        Assert.AreEqual(t3[14].pointarray[2], new Point64(370,320,0));
        Assert.AreEqual(t3[14].pointarray[3], new Point64(370,260,0));
        Assert.AreEqual(t3[14].pointarray[4], new Point64(430,260,0));
        Assert.AreEqual(t3[14].pointarray[5], new Point64(430,220,0));
        Assert.AreEqual(t3[14].pointarray[6], new Point64(330,220,0));
        Assert.AreEqual(t3[15].pointarray.Count, 7);
        Assert.AreEqual(t3[15].pointarray[0], new Point64(330,330,0));
        Assert.AreEqual(t3[15].pointarray[1], new Point64(330,430,0));
        Assert.AreEqual(t3[15].pointarray[2], new Point64(370,430,0));
        Assert.AreEqual(t3[15].pointarray[3], new Point64(370,370,0));
        Assert.AreEqual(t3[15].pointarray[4], new Point64(430,370,0));
        Assert.AreEqual(t3[15].pointarray[5], new Point64(430,330,0));
        Assert.AreEqual(t3[15].pointarray[6], new Point64(330,330,0));
    }

    [Test]
    public static void consistency_from_oasis()
    {
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(baseDir + "/consistency/triage.oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);

        string outFile = baseDir + "/out/c3_consistency_from_oas.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcOAS, outFile);
        gw.save();
        Assert.True(File.Exists(outFile));

        string outFile2 = baseDir + "/out/c3_consistency_from_oas.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcOAS, outFile2);
        ow.save();
        Assert.True(File.Exists(outFile2));
    }

    [Test]
    public static void consistency_from_gds()
    {
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(baseDir + "/consistency/c3_consistency.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);

        string outFile = baseDir + "/out/c3_consistency_from_gds.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcGDS, outFile);
        gw.save();
        Assert.True(File.Exists(outFile));

        string outFile2 = baseDir + "/out/c3_consistency_from_gds.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcGDS, outFile2);
        ow.save();
        Assert.True(File.Exists(outFile2));
    }

    [Test]
    public static void test_cell_export_complex()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested_2.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");
        Assert.AreEqual(gcGDS.activeStructure, 0);

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");
        Assert.AreEqual(gcGDS.activeLD, 0);

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.AreEqual(t1.Count, 1);
        Assert.AreEqual(t1[0].layer_nr, 1);
        Assert.AreEqual(t1[0].datatype_nr, 0);
        Assert.AreEqual(t1[0].pointarray.Count, 7);
        Assert.AreEqual(t1[0].pointarray[0], new Point64(0,0,0));
        Assert.AreEqual(t1[0].pointarray[1], new Point64(0,100,0));
        Assert.AreEqual(t1[0].pointarray[2], new Point64(40,100,0));
        Assert.AreEqual(t1[0].pointarray[3], new Point64(40,40,0));
        Assert.AreEqual(t1[0].pointarray[4], new Point64(100,40,0));
        Assert.AreEqual(t1[0].pointarray[5], new Point64(100,0,0));
        Assert.AreEqual(t1[0].pointarray[6], new Point64(0,0,0));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.AreEqual(gcGDS.activeStructure, 1);

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.AreEqual(t2.Count, 6);
        Assert.AreEqual(t2[0].layer_nr, 1);
        Assert.AreEqual(t2[0].datatype_nr, 0);
        Assert.AreEqual(t2[0].pointarray.Count, 7);
        Assert.AreEqual(t2[0].pointarray[0], new Point64(0,0,0));
        Assert.AreEqual(t2[0].pointarray[1], new Point64(0,100,0));
        Assert.AreEqual(t2[0].pointarray[2], new Point64(40,100,0));
        Assert.AreEqual(t2[0].pointarray[3], new Point64(40,40,0));
        Assert.AreEqual(t2[0].pointarray[4], new Point64(100,40,0));
        Assert.AreEqual(t2[0].pointarray[5], new Point64(100,0,0));
        Assert.AreEqual(t2[0].pointarray[6], new Point64(0,0,0));
        Assert.AreEqual(t2[1].layer_nr, 1);
        Assert.AreEqual(t2[1].datatype_nr, 0);
        Assert.AreEqual(t2[1].pointarray.Count, 7);
        Assert.AreEqual(t2[1].pointarray[0], new Point64(0,110,0));
        Assert.AreEqual(t2[1].pointarray[1], new Point64(0,210,0));
        Assert.AreEqual(t2[1].pointarray[2], new Point64(40,210,0));
        Assert.AreEqual(t2[1].pointarray[3], new Point64(40,150,0));
        Assert.AreEqual(t2[1].pointarray[4], new Point64(100,150,0));
        Assert.AreEqual(t2[1].pointarray[5], new Point64(100,110,0));
        Assert.AreEqual(t2[1].pointarray[6], new Point64(0,110,0));
        Assert.AreEqual(t2[2].layer_nr, 1);
        Assert.AreEqual(t2[2].datatype_nr, 0);
        Assert.AreEqual(t2[2].pointarray.Count, 7);
        Assert.AreEqual(t2[2].pointarray[0], new Point64(110,0,0));
        Assert.AreEqual(t2[2].pointarray[1], new Point64(110,100,0));
        Assert.AreEqual(t2[2].pointarray[2], new Point64(150,100,0));
        Assert.AreEqual(t2[2].pointarray[3], new Point64(150,40,0));
        Assert.AreEqual(t2[2].pointarray[4], new Point64(210,40,0));
        Assert.AreEqual(t2[2].pointarray[5], new Point64(210,0,0));
        Assert.AreEqual(t2[2].pointarray[6], new Point64(110,0,0));
        Assert.AreEqual(t2[3].layer_nr, 1);
        Assert.AreEqual(t2[3].datatype_nr, 0);
        Assert.AreEqual(t2[3].pointarray.Count, 7);
        Assert.AreEqual(t2[3].pointarray[0], new Point64(110,110,0));
        Assert.AreEqual(t2[3].pointarray[1], new Point64(110,210,0));
        Assert.AreEqual(t2[3].pointarray[2], new Point64(150,210,0));
        Assert.AreEqual(t2[3].pointarray[3], new Point64(150,150,0));
        Assert.AreEqual(t2[3].pointarray[4], new Point64(210,150,0));
        Assert.AreEqual(t2[3].pointarray[5], new Point64(210,110,0));
        Assert.AreEqual(t2[3].pointarray[6], new Point64(110,110,0));
        Assert.AreEqual(t2[4].layer_nr, 1);
        Assert.AreEqual(t2[4].datatype_nr, 0);
        Assert.AreEqual(t2[4].pointarray.Count, 5);
        Assert.AreEqual(t2[4].pointarray[0], new Point64(60,173,0));
        Assert.AreEqual(t2[4].pointarray[1], new Point64(60,201,0));
        Assert.AreEqual(t2[4].pointarray[2], new Point64(83,201,0));
        Assert.AreEqual(t2[4].pointarray[3], new Point64(83,173,0));
        Assert.AreEqual(t2[4].pointarray[4], new Point64(60,173,0));
        Assert.AreEqual(t2[5].layer_nr, 2);
        Assert.AreEqual(t2[5].datatype_nr, 0);
        Assert.AreEqual(t2[5].pointarray.Count, 5);
        Assert.AreEqual(t2[5].pointarray[0], new Point64(52,67,0));
        Assert.AreEqual(t2[5].pointarray[1], new Point64(52,93,0));
        Assert.AreEqual(t2[5].pointarray[2], new Point64(103,93,0));
        Assert.AreEqual(t2[5].pointarray[3], new Point64(103,67,0));
        Assert.AreEqual(t2[5].pointarray[4], new Point64(52,67,0));

        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        List<GCPolygon> t3a = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(layer: 1, datatype: 0);
        List<GCPolygon> t3b = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(layer: 2, datatype: 0);

        List<GCPolygon> t4 = gcGDS.convertToPolygons();
        List<GCPolygon> t4a = gcGDS.convertToPolygons(activeLDOnly: true);
        List<GCPolygon> t4b = gcGDS.convertToPolygons(layer: 2, datatype: 0);

        // Ensure we have equivalence per expectations.
        Assert.AreEqual(Utils.GetSHA256Hash(t3), Utils.GetSHA256Hash(t4));
        Assert.AreNotEqual(Utils.GetSHA256Hash(t3), Utils.GetSHA256Hash(t3a));
        Assert.AreNotEqual(Utils.GetSHA256Hash(t3), Utils.GetSHA256Hash(t3b));
        Assert.AreEqual(Utils.GetSHA256Hash(t3a), Utils.GetSHA256Hash(t4a));
        Assert.AreEqual(Utils.GetSHA256Hash(t3b), Utils.GetSHA256Hash(t4b));
    }

    [Test]
    public static void test_cellrefarray_basic()
    {
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // The array is in cell 'a'
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        PathsD geo = gcGDS.points(flatten: true);
        Assert.AreEqual(geo.Count, 4);
        
        // We have floats, so this gets a little more awkward.
        Assert.LessOrEqual(Math.Abs(geo[0][0].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][0].y - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][1].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][1].y - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][2].x - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][2].y - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][3].x - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][3].y - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][4].x - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][4].y - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][5].x - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][5].y - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][6].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[0][6].y - 0), 1E-13);

        int y_adjust = 110;
        Assert.LessOrEqual(Math.Abs(geo[1][0].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][0].y - (0 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][1].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][1].y - (100 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][2].x - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][2].y - (100 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][3].x - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][3].y - (40 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][4].x - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][4].y - (40 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][5].x - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][5].y - 0 - y_adjust), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][6].x - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[1][6].y - (0 + y_adjust)), 1E-13);

        int x_adjust = 110;
        Assert.LessOrEqual(Math.Abs(geo[2][0].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][0].y - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][1].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][1].y - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][2].x - (40 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][2].y - 100), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][3].x - (40 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][3].y - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][4].x - (100 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][4].y - 40), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][5].x - (100 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][5].y - 0), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][6].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[2][6].y - 0), 1E-13);

        Assert.LessOrEqual(Math.Abs(geo[3][0].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][1].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][2].x - (40 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][3].x - (40 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][4].x - (100 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][5].x - (100 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][6].x - (0 + x_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][0].y - (0 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][1].y - (100 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][2].y - (100 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][3].y - (40 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][4].y - (40 + y_adjust)), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][5].y - 0 - y_adjust), 1E-13);
        Assert.LessOrEqual(Math.Abs(geo[3][6].y - (0 + y_adjust)), 1E-13);
    }

    [Test]
    public static void test_cellrefarray_nested()
    {
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested_2.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();

        // The array is in cell 'a'
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        PathsD geo2 = gcGDS.points(flatten: true);
        // Use previously computed hash to check that our long list is aligned with expectations.
        Assert.AreEqual(Utils.GetSHA256Hash(geo2), "9Y7vZPon4aFRZ07yD8a3dzuaTq9n1mMnHJaE89C3fb0=");
    }
    
    [Test]
    public static void defineAndWrite_Polygon()
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
        Path64 poly = Helper.initedPath64(6);
        poly[0] = new (0, 0);
        poly[1] = new (0, 20);
        poly[2] = new (10, 20);
        poly[3] = new (10, 10);
        poly[4] = new (20, 10);
        poly[5] = new (20, 0);

        gcell.addPolygon(poly, 1, 0);

        // triangle
        poly = Helper.initedPath64(3);
        poly[0] = new (0, 0);
        poly[1] = new (10, 20);
        poly[2] = new (20, 0);

        gcell.addPolygon(poly, 2, 0);

        // pentagram
        poly = Helper.initedPath64(5);
        poly[0] = new (5, 0);
        poly[1] = new (0, 10);
        poly[2] = new (10, 20);
        poly[3] = new (20, 10);
        poly[4] = new (15, 0);

        gcell.addPolygon(poly, 3, 0);

        // trapezoid
        poly = Helper.initedPath64(4);
        poly[0] = new (0, 0);
        poly[1] = new (5, 20);
        poly[2] = new (15, 20);
        poly[3] = new (20, 0);

        gcell.addPolygon(poly, 4, 0);

        // parallelogram
        poly = Helper.initedPath64(4);
        poly[0] = new (0, 0);
        poly[1] = new (10, 20);
        poly[2] = new (20, 20);
        poly[3] = new (10, 0);

        gcell.addPolygon(poly, 5, 0);

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_polygon.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();

        Assert.True(File.Exists(gdsFile));
        
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        // L
        Assert.AreEqual(true, cell_gds.elementList[0].isPolygon());
        Assert.AreEqual(1, cell_gds.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[0].datatype_nr);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(7, polys_gds[0].pointarray.Count);

        // Triangle
        Assert.AreEqual(true, cell_gds.elementList[1].isPolygon());
        Assert.AreEqual(2, cell_gds.elementList[1].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[1].datatype_nr);
        polys_gds = cell_gds.elementList[1].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(4, polys_gds[0].pointarray.Count);

        // Pentagram
        Assert.AreEqual(true, cell_gds.elementList[2].isPolygon());
        Assert.AreEqual(3, cell_gds.elementList[2].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[2].datatype_nr);
        polys_gds = cell_gds.elementList[2].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(6, polys_gds[0].pointarray.Count);

        // Trapezoid
        Assert.AreEqual(true, cell_gds.elementList[3].isPolygon());
        Assert.AreEqual(4, cell_gds.elementList[3].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[3].datatype_nr);
        polys_gds = cell_gds.elementList[3].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(5, polys_gds[0].pointarray.Count);

        // Parallelogram
        Assert.AreEqual(true, cell_gds.elementList[4].isPolygon());
        Assert.AreEqual(5, cell_gds.elementList[4].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[4].datatype_nr);
        polys_gds = cell_gds.elementList[4].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(5, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_polygon.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        // L
        Assert.AreEqual(true, cell_oas.elementList[0].isPolygon());
        Assert.AreEqual(1, cell_oas.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[0].datatype_nr);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(7, polys_oas[0].pointarray.Count);

        // Triangle
        Assert.AreEqual(true, cell_oas.elementList[1].isPolygon());
        Assert.AreEqual(2, cell_oas.elementList[1].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[1].datatype_nr);
        polys_oas = cell_oas.elementList[1].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(4, polys_oas[0].pointarray.Count);

        // Pentagram
        Assert.AreEqual(true, cell_oas.elementList[2].isPolygon());
        Assert.AreEqual(3, cell_oas.elementList[2].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[2].datatype_nr);
        polys_oas = cell_oas.elementList[2].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(6, polys_oas[0].pointarray.Count);

        // Trapezoid
        Assert.AreEqual(true, cell_oas.elementList[3].isPolygon());
        Assert.AreEqual(4, cell_oas.elementList[3].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[3].datatype_nr);
        polys_oas = cell_oas.elementList[3].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(5, polys_oas[0].pointarray.Count);

        // Parallelogram
        Assert.AreEqual(true, cell_oas.elementList[4].isPolygon());
        Assert.AreEqual(5, cell_oas.elementList[4].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[4].datatype_nr);
        polys_oas = cell_oas.elementList[4].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(5, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_Path()
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

        Path64 path = Helper.initedPath64(5);
        path[0] = new (0, 0);
        path[1] = new (0, 10);
        path[2] = new (20, 10);
        path[3] = new (20, 40);
        path[4] = new (0, 40);

        gcell.addPath(path, 1, 0);
        gcell.elementList[^1].setWidth(5); // note that Oasis only supports factors of 2, so this gets rounded down to make a 4 unit path at the writer (for Oasis).
        gcell.elementList[^1].setCap(2); // extend the path by half the width at the line ends.

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_path.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.AreEqual(true, cell_gds.elementList[0].isPath());
        Assert.AreEqual(1, cell_gds.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[0].datatype_nr);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(12, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_path.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.AreEqual(true, cell_oas.elementList[0].isPath());
        Assert.AreEqual(1, cell_oas.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[0].datatype_nr);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(12, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_Circle()
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

        gcell.addCircle(1, 0, new (10, 10), 5.0);

        g.setDrawing(drawing_);
        g.setValid(true);

        // GDS looks janky for circles; no special category unlike the Oasis format. Snaps to grid and looks ugly.
        string gdsFile = outDir + "simple_circle.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.AreEqual(true, cell_gds.elementList[0].isPolygon());
        Assert.AreEqual(1, cell_gds.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_gds.elementList[0].datatype_nr);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(3600, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_circle.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.AreEqual(true, cell_oas.elementList[0].isPolygon());
        Assert.AreEqual(1, cell_oas.elementList[0].layer_nr);
        Assert.AreEqual(0, cell_oas.elementList[0].datatype_nr);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(3600, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_Cellref()
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

        Path64 poly = Helper.initedPath64(6);
        poly[0] = new (0, 0);
        poly[1] = new (0, 20);
        poly[2] = new (10, 20);
        poly[3] = new (10, 10);
        poly[4] = new (20, 10);
        poly[5] = new (20, 0);

        gcell.addPolygon(poly, 1, 0);

        bool mirror_x = false;
        gcell = drawing_.addCell();
        gcell.addCellref();

        gcell.elementList[^1].setPos(new (10, 0));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test"));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[^1].setMirrorx();
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

        poly = Helper.initedPath64(6);
        poly[0] = new (0, 0);
        poly[1] = new (0, 30);
        poly[2] = new (10, 30);
        poly[3] = new (10, 10);
        poly[4] = new (20, 10);
        poly[5] = new (20, 0);

        gcell.addPolygon(poly, 1, 0);

        mirror_x = true;
        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new (20, 20));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test2"));
        gcell.elementList[^1].setName("test2");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[^1].setMirrorx();
                break;
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellref.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellref");

        Assert.AreEqual(true, cell_gds.elementList[^1].isCellref());
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(7, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_cellref.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellref");

        Assert.AreEqual(true, cell_oas.elementList[^1].isCellref());
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(7, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_CellrefArray()
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

        Path64 poly = Helper.initedPath64(6);
        poly[0] = new (0, 0);
        poly[1] = new (0, 20);
        poly[2] = new (10, 20);
        poly[3] = new (10, 10);
        poly[4] = new (20, 10);
        poly[5] = new (20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new (0, 0);
        array[1] = new (100, 0);
        array[2] = new (0, 80);

        bool mirror_x = false;
        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (0, 0));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);
        switch (mirror_x)
        {
            case true:
                gcell.elementList[^1].setMirrorx();
                break;
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray");

        Assert.AreEqual(true, cell_gds.elementList[^1].isCellrefArray());
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.AreEqual(16, polys_gds.Count);
        Assert.AreEqual(7, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_cellrefarray.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray");

        Assert.AreEqual(true, cell_oas.elementList[^1].isCellrefArray());
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.AreEqual(16, polys_oas.Count);
        Assert.AreEqual(7, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_Text()
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

        gcell.addText(1, 0, new (10, 10), "Text");

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_text.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.AreEqual(true, cell_gds.elementList[^1].isText());

        string oasFile = outDir + "simple_text.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.AreEqual(true, cell_oas.elementList[^1].isText());
    }

    [Test]
    public static void defineAndWrite_CTrapezoid()
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
            Path64 pa = Helper.initedPath64(5);

            int[,] coords = i switch
            {
                0 => new int[4, 4]
                {
                    {0, 0, 0, 0}, // x=0*w+0*h, y=0*w+0*h ...
                    {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}
                },
                1 => new[,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                2 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                3 => new[,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                4 => new[,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
                5 => new[,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                6 => new [,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {1, 0, 0, 1}, {1, -1, 0, 0}},
                7 => new [,] {{0, 1, 0, 0}, {0, 0, 0, 1}, {1, -1, 0, 1}, {1, 0, 0, 0}},
                8 => new [,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
                9 => new [,] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                10 => new [,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
                11 => new [,] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                12 => new [,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 1, 0}},
                13 => new [,] {{0, 0, 1, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                14 => new [,] {{0, 0, 0, 0}, {0, 0, -1, 1}, {1, 0, 0, 1}, {1, 0, 1, 0}},
                15 => new [,] {{0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, -1, 1}, {1, 0, 0, 0}},
                16 => new [,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
                17 => new [,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
                18 => new [,] {{0, 0, 0, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 0, 0}},
                19 => new [,] {{0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}},
                20 => new [,] {{0, 0, 0, 0}, {0, 1, 0, 1}, {0, 2, 0, 0}, {0, 0, 0, 0}},
                21 => new [,] {{0, 0, 0, 1}, {0, 2, 0, 1}, {0, 1, 0, 0}, {0, 0, 0, 1}},
                22 => new [,] {{0, 0, 0, 0}, {0, 0, 2, 0}, {1, 0, 1, 0}, {0, 0, 0, 0}},
                23 => new [,] {{1, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 2, 0}, {1, 0, 0, 0}},
                24 => new [,] {{0, 0, 0, 0}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}},
                25 => new [,] {{0, 0, 0, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 0, 0}},
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

                pa[pt] = new (x, y);

                if (x > w)
                {
                    w = x;
                }
                if (y > h)
                {
                    h = y;
                }
            }

            pa[^1] = new (pa[0]);

            gcell.addPolygon(pa, i + 1, 0);
        }


        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_ctrapezoid.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.AreEqual(true, cell_gds.elementList[^1].isPolygon());
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(5, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_ctrapezoid.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.AreEqual(true, cell_oas.elementList[^1].isPolygon());
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(5, polys_oas[0].pointarray.Count);
    }

    [Test]
    public static void defineAndWrite_Trapezoid()
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

        Path64 pa = Helper.initedPath64(5);

        switch (trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new (x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new (x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new (x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new (x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new (x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new (x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new (x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new (x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new (pa[0]);

        gcell.addPolygon(pa, 1, 0);

        trapezoid_orientation = true;

        switch (trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new (x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new (x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new (x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new (x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new (x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new (x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new (x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new (x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new (pa[0]);

        gcell.addPolygon(pa, 2, 0);

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_trapezoid.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.AreEqual(true, cell_gds.elementList[^1].isPolygon());
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_gds.Count);
        Assert.AreEqual(5, polys_gds[0].pointarray.Count);

        string oasFile = outDir + "simple_trapezoid.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.AreEqual(true, cell_oas.elementList[^1].isPolygon());
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.AreEqual(1, polys_oas.Count);
        Assert.AreEqual(5, polys_oas[0].pointarray.Count);
    }

    private static void defineAndWrite_Box_LayerDataType()
    {
        int[] values = { 1, 128, 255, 256, 511, 512, 16383, 16384, 32767, 32768 };
        foreach (int layer in values)
        {
            int datatype = 0;
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

            string gdsFile = outDir + "L" + layer + "D" + datatype + "_box.gds";
            if (File.Exists(gdsFile))
            {
                File.Delete(gdsFile);
            }
            gds.gdsWriter gw = new(g, gdsFile);
            gw.save();
            Assert.True(File.Exists(gdsFile));

            string oasFile = outDir + "L" + layer + "D" + datatype + "_box.oas";
            if (File.Exists(oasFile))
            {
                File.Delete(oasFile);
            }
            oasis.oasWriter ow = new(g, oasFile);
            ow.save();
            Assert.True(File.Exists(oasFile));
        }
        
        foreach (int datatype in values)
        {
            int layer = 1;
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

            string gdsFile = outDir + "L" + layer + "D" + datatype + "_box.gds";
            if (File.Exists(gdsFile))
            {
                File.Delete(gdsFile);
            }
            gds.gdsWriter gw = new(g, gdsFile);
            gw.save();
            Assert.True(File.Exists(gdsFile));

            string oasFile = outDir + "L" + layer + "D" + datatype + "_box.oas";
            if (File.Exists(oasFile))
            {
                File.Delete(oasFile);
            }
            oasis.oasWriter ow = new(g, oasFile);
            ow.save();
            Assert.True(File.Exists(oasFile));
        }
    }

    [Test]
    public static void defineAndWrite_manyCellrefs()
    {
        int edge = 4;
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

            Path64 poly = Helper.initedPath64(6);
            poly[0] = new (0, 0);
            poly[1] = new (0, 20);
            poly[2] = new (10, 20);
            poly[3] = new (10, 10);
            poly[4] = new (20, 10);
            poly[5] = new (20, 0);

            gcell.addPolygon(poly, 1, 0);

            master_gcell.addCellref();
            GCElement cellref = master_gcell.elementList[^1];
            cellref.setPos(new (40 * (i % edge), 40 * Math.Floor((double)i / edge)));
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

        string gdsFile = outDir + edge + "_cellref.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.True(File.Exists(gdsFile));

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.AreEqual(gcGDS.isValid(), true);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        for (int i = 0; i < edge * edge; i++)
        {
            GCCell cell_gds = drawing_gds.findCell("test" + i);

            Assert.AreEqual(true, cell_gds.elementList[^1].isPolygon());
            List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
            Assert.AreEqual(1, polys_gds.Count);
            Assert.AreEqual(7, polys_gds[0].pointarray.Count);
        }

        string oasFile = outDir + edge + "_cellref.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.True(File.Exists(oasFile));
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.AreEqual(gcOAS.isValid(), true);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        for (int i = 0; i < edge * edge; i++)
        {
            GCCell cell_oas = drawing_oas.findCell("test" + i);

            Assert.AreEqual(true, cell_oas.elementList[^1].isPolygon());
            List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
            Assert.AreEqual(1, polys_oas.Count);
            Assert.AreEqual(7, polys_oas[0].pointarray.Count);
        }
    }
}