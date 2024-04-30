using System.Buffers.Binary;
using Clipper2Lib;
using geoCoreLib;
using utility;

namespace UnitTests;

public class GeoCoreTests
{
    static string baseDir = "/d/development/DesignLibs_GPL/geocore_test/";
    static string outDir = "/d/development/geocore_out/";

    // [SetUp]
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
        defineAndWrite_CellrefArray1();
        defineAndWrite_CellrefArray2();
        defineAndWrite_CellrefArray3();
        defineAndWrite_CellrefArray4();
        defineAndWrite_CellrefArray5();
        defineAndWrite_CellrefArray1mx();
        defineAndWrite_CellrefArray2mx();
        defineAndWrite_CellrefArray3mx();
        defineAndWrite_CellrefArray4mx();
        defineAndWrite_CellrefArray5mx();
        defineAndWrite_CellrefArray_irregular();
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
    public static void f_rep_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isPolygon(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(1200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(1100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(1100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(1200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(1200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(1100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(1100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(-2000));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(-2000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(-1800));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(-1800));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(-1900));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(-1900));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(-1800));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(-1800));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(-1900));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(-1900));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(-2000));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        // OASIS file appears to be invalid per KLayout; bug raised with gdstk
        /*
        string oasFile = baseDir + "gdstk_reference/f_rep.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        */
    }
    
    [Test]
    public static void f_rep2_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep2.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));
        
        int polyIndex = 0;
        Assert.That(cell_gds.elementList[polyIndex].isPolygon(), Is.True);
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        Assert.That(cell_gds.elementList[polyIndex].isPolygon(), Is.True);
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(4000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(4000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(3800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(3800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(3600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(3600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(3200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(3400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(3400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(3100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(3000));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(3000));
        
        string oasFile = baseDir + "gdstk_reference/f_rep2.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        Assert.That(cell_oas.elementList[polyIndex].isPolygon(), Is.True);
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(3000));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(3000));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(3000));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(4000));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(3200));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(4000));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(3200));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(3800));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(3100));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(3800));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(3100));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(3600));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(3200));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(3600));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(3200));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(3400));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(3100));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(3400));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(3100));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(3000));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(3000));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(3000));
    }
    
    [Test]
    public static void f_rep3_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep3.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons();
        
        int polyIndex = 0;
        int x_pitch = 5000;
        int y_pitch = 5000;
        for (int col = 0; col < 2; col++)
        {
            int x_offset = col * x_pitch;
            for (int row = 0; row < 3; row++)
            {
                int y_offset = row * y_pitch;
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
                polyIndex++;
            }
        }

        string oasFile = baseDir + "gdstk_reference/f_rep3.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons();

        polyIndex = 0;
        for (int row = 0; row < 3; row++)
        {
            int y_offset = row * y_pitch;
            for (int col = 0; col < 2; col++)
            {
                int x_offset = col * x_pitch;
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
                polyIndex++;
            }
        }
    }
    
    [Test]
    public static void f_rep4_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep4.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons();
        
        int polyIndex = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        int x_offset = 500;
        int y_offset = 1000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        string oasFile = baseDir + "gdstk_reference/f_rep4.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons();
        
        polyIndex = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        x_offset = 500;
        y_offset = 1000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void f_rep5_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep5.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons();
        
        int polyIndex = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0));

        polyIndex++;
        int x_offset = 4000;
        int y_offset = 4000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        string oasFile = baseDir + "gdstk_reference/f_rep5.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons();

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 4000;
        y_offset = 4000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void rect_rep_test()
    {
        string gdsFile = baseDir + "gdstk_reference/rectangle_rep.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000));

        polyIndex++;
        int x_offset = 1000;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        string oasFile = baseDir + "gdstk_reference/rectangle_rep.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcGDS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_oas.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 1000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void rect_rep2_test()
    {
        string gdsFile = baseDir + "gdstk_reference/rectangle_rep2.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_gds.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));

        int polyIndex = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000));

        polyIndex++;
        int x_offset = 3000;
        int y_offset = 3000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        string oasFile = baseDir + "gdstk_reference/rectangle_rep2.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcGDS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_oas.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void rect_rep3_test()
    {
        string gdsFile = baseDir + "gdstk_reference/rectangle_rep3.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(6));
        for (int i = 0; i < 6; i++)
        {
            Assert.That(cell_gds.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(6));

        int polyIndex = 0;
        int x_pitch = 5000;
        int y_pitch = 5000;
        for (int col = 0; col < 2; col++)
        {
            int x_offset = col * x_pitch;
            for (int row = 0; row < 3; row++)
            {
                int y_offset = row * y_pitch;
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
                polyIndex++;
            }
        }
        
        string oasFile = baseDir + "gdstk_reference/rectangle_rep3.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcGDS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(6));
        for (int i = 0; i < 6; i++)
        {
            Assert.That(cell_oas.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(6));

        polyIndex = 0;
        for (int col = 0; col < 2; col++)
        {
            int x_offset = col * x_pitch;
            for (int row = 0; row < 3; row++)
            {
                int y_offset = row * y_pitch;
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
                polyIndex++;
            }
        }
    }
    
    [Test]
    public static void rect_rep4_test()
    {
        string gdsFile = baseDir + "gdstk_reference/rectangle_rep4.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 500;
        y_offset = 1000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        string oasFile = baseDir + "gdstk_reference/rectangle_rep4.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_oas.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 500;
        y_offset = 1000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void rect_rep5_test()
    {
        string gdsFile = baseDir + "gdstk_reference/rectangle_rep5.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Base");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_gds.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 4000;
        y_offset = 4000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        string oasFile = baseDir + "gdstk_reference/rectangle_rep5.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_oas.elementList[i].isBox(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 4000;
        y_offset = 4000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void ref_f_rep_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        string oasFile = baseDir + "gdstk_reference/ref_f_rep.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void ref_f_rep2_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep2.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        // This one gets handled as a 'regular' repetition.
        // This internally then gets mapped to a cellrefarray
        string oasFile = baseDir + "gdstk_reference/ref_f_rep2.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
        
        // This one gets handled as an 'explicit' repetition. The file was simply passed through KLayout and ended up
        // with this different representation.
        oasFile = baseDir + "gdstk_reference/ref_f_rep2_klayout.oas";
        GeoCoreHandler gH_OAS_kl = new();
        gH_OAS_kl.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS_kl = gH_OAS_kl.getGeo();
        Assert.That(gcOAS_kl.isValid(), Is.True);
        
        GCDrawingfield drawing_oas_kl = gcOAS_kl.getDrawing();
        GCCell cell_oas_kl = drawing_oas_kl.findCell("Ref");
        Assert.That(cell_oas_kl.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_oas_kl.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_oas_kl = cell_oas_kl.convertToPolygons();
        Assert.That(polys_oas_kl.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas_kl[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas_kl[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_oas_kl[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas_kl[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas_kl[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void ref_f_rep3_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep3.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(6));

        int polyIndex = 0;
        int x_pitch = 5000;
        int y_pitch = 5000;
        for (int row = 0; row < 3; row++)
        {
            int y_offset = row * y_pitch;
            for (int col = 0; col < 2; col++)
            {
                int x_offset = col * x_pitch;
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
                polyIndex++;
            }
        }

        string oasFile = baseDir + "gdstk_reference/ref_f_rep3.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(6));

        polyIndex = 0;
        for (int row = 0; row < 3; row++)
        {
            int y_offset = row * y_pitch;
            for (int col = 0; col < 2; col++)
            {
                int x_offset = col * x_pitch;
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
                polyIndex++;
            }
        }
    }
    
    [Test]
    public static void ref_f_rep4_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep4.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 500;
        y_offset = 1000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        string oasFile = baseDir + "gdstk_reference/ref_f_rep4.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 500;
        y_offset = 1000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 1500;
        y_offset = 500;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void ref_f_rep5_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep5.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 4000;
        y_offset = 4000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
        
        // The OASIS reader captures the cellref array as a more efficient representation.
        string oasFile = baseDir + "gdstk_reference/ref_f_rep5.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));

        polyIndex++;
        x_offset = 4000;
        y_offset = 4000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(800 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(600 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].X, Is.EqualTo(200 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[7].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[8].Y, Is.EqualTo(400 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].X, Is.EqualTo(100 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[9].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[10].Y, Is.EqualTo(0 + y_offset));
    }
    
    [Test]
    public static void ref_rect_rep_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(4));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 1000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        x_offset = -2000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 1000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void ref_rect_rep2_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep2.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        for (int i = 0; i < 2; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(2));

        int polyIndex = 0;
        int x_offset = 0;
        int y_offset = 0;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
        
        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep2.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(2));

        polyIndex = 0;
        x_offset = 0;
        y_offset = 0;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));

        polyIndex++;
        x_offset = 3000;
        y_offset = 3000;
        Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
        Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
    }
    
    [Test]
    public static void ref_rect_rep3_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep3.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_gds = cell_gds.convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(6));

        int polyIndex = 0;
        int x_pitch = 5000;
        int y_pitch = 5000;
        for (int row = 0; row < 3; row++)
        {
            int y_offset = row * y_pitch;
            for (int col = 0; col < 2; col++)
            {
                int x_offset = col * x_pitch;
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(5));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
                polyIndex++;
            }
        }

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep3.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        
        GCDrawingfield drawing_oas = gcGDS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Ref");
        Assert.That(cell_oas.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_oas.elementList[i].isCellrefArray(), Is.True);
        }
        List <GCPolygon> polys_oas = cell_oas.convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(6));

        polyIndex = 0;
        for (int row = 0; row < 3; row++)
        {
            int y_offset = row * y_pitch;
            for (int col = 0; col < 2; col++)
            {
                int x_offset = col * x_pitch;
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(5));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(1000 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(1000 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(0 + y_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + x_offset));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(1000 + y_offset));
                polyIndex++;
            }
        }
    }
    
    [Test]
    public static void ref_rect_rep4_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep4.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep4.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
    }
    
    [Test]
    public static void ref_rect_rep5_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep5.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep5.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
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
            Assert.That(e.Message, Is.EqualTo("Provided GeoCore instance is not marked as valid"));
        }
        try
        {
            oasis.oasWriter ow = new(g, outDir + "/should_fail.oas");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.That(e.Message, Is.EqualTo("Provided GeoCore instance is not marked as valid"));
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
            Assert.That(e.Message, Is.EqualTo("Provided GeoCore instance is not marked as valid"));
        }
        try
        {
            oasis.oasWriter ow = new(g, outDir + "/should_fail.oas");
            throw new("Failed to fail");
        }
        catch (Exception e)
        {
            Assert.That(e.Message, Is.EqualTo("Provided GeoCore instance is not marked as valid"));
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
                Assert.That(drawing_.accyear, Is.EqualTo(2018));

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

                Assert.That(gcell.elementList[^1].isBox(), Is.True);

                g.setDrawing(drawing_);
                g.setValid(true);

                if (File.Exists(outDir + "/" + filename + ".gds"))
                {
                    File.Delete(outDir + "/" + filename + ".gds");
                }

                gds.gdsWriter gw = new(g, outDir + "/" + filename + ".gds");
                gw.save();
                Assert.That(File.Exists(outDir + "/" + filename + ".gds"), Is.True);

                if (File.Exists(outDir + "/" + filename + ".oas"))
                {
                    File.Delete(outDir + "/" + filename + ".oas");
                }

                oasis.oasWriter ow = new(g, outDir + "/" + filename + ".oas");
                ow.save();
                Assert.That(File.Exists(outDir + "/" + filename + ".oas"), Is.True);

                // Load the files in to see what we have.

                GeoCoreHandler gH_GDS = new();
                gH_GDS.updateGeoCoreHandler(outDir + "/" + filename + ".gds", GeoCore.fileType.gds);
                GeoCore gcGDS = gH_GDS.getGeo();
                Assert.That(gcGDS.isValid(), Is.True);

                GCDrawingfield drawing_gds = gcGDS.getDrawing();
                drawing_gds.databaseunits = 1000 * scale;
                drawing_gds.userunits = 0.001 / scale;
                GCCell cell_gds = drawing_gds.findCell("test");
                int elementCount = cell_gds.elementList.Count;
                Assert.That(elementCount, Is.EqualTo(1));
                // Assert.True(cell_gds.elementList[^1].isBox());

                string out_filename = outDir + "/" + filename + "_resave_from_gds.gds";
                save_gdsii(gcGDS, out_filename);

                out_filename = outDir + "/" + filename + "_resave_from_gds.oas";
                save_oasis(gcGDS, out_filename);

                GeoCoreHandler gH_OAS = new();
                gH_OAS.updateGeoCoreHandler(outDir + "/" + filename + ".oas", GeoCore.fileType.oasis);
                GeoCore gcOAS = gH_OAS.getGeo();
                Assert.That(gcOAS.isValid(), Is.True);

                GCDrawingfield drawing_oas = gcOAS.getDrawing();
                drawing_oas.databaseunits = 1000 * scale;
                drawing_oas.userunits = 0.001 / scale;
                GCCell cell_oas = drawing_oas.findCell("test");
                elementCount = cell_oas.elementList.Count;
                Assert.That(elementCount, Is.EqualTo(1));
                Assert.That(cell_oas.elementList[^1].isBox(), Is.True);

                out_filename = outDir + "/" + filename + "_resave_from_oas.gds";
                save_gdsii(gcOAS, out_filename);

                out_filename = outDir + "/" + filename + "_resave_from_oas.oas";
                save_oasis(gcOAS, out_filename);
            }
        }
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
        Assert.That(drawing_.accyear, Is.EqualTo(2018));

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
        Assert.That(drawing_.findCell("test"), Is.EqualTo(gcell));
       
        g.setDrawing(drawing_);
        g.setValid(true);
        
        if (File.Exists(outDir + "poly10_circle11.gds"))
        {
            File.Delete(outDir + "poly10_circle11.gds");
        }
        gds.gdsWriter gw = new(g, outDir + "poly10_circle11.gds");
        gw.save();
        Assert.That(File.Exists(outDir + "poly10_circle11.gds"), Is.True);

        if (File.Exists(outDir + "poly10_circle11.oas"))
        {
            File.Delete(outDir + "poly10_circle11.oas");
        }
        oasis.oasWriter ow = new(g, outDir + "poly10_circle11.oas");
        ow.save();
        Assert.That(File.Exists(outDir + "poly10_circle11.oas"), Is.True);
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
        Assert.That(drawing_.accyear, Is.EqualTo(2018));

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

        PathD poly1 = Clipper.MakePath(new double[]
        {
            60, 460,
            79, 460,
            79, 470,
            60, 470,
            60, 460
        });
        PathD poly2 = Clipper.MakePath(new double[]
        {
            70, 459,
            79, 459,
            79, 460,
            70, 460,
            70, 459,
        });
        PathD poly3 = Clipper.MakePath(new double[]
        {
            70, 450,
            80, 450,
            80, 459,
            70, 459,
            70, 450
        });
        
        gcell.addPolygon(Clipper.ScalePath64(poly1, 1), 1, 0);
        gcell.addPolygon(Clipper.ScalePath64(poly2, 1), 1, 0);
        gcell.addPolygon(Clipper.ScalePath64(poly3, 1), 1, 0);

        // gcell.addPolygon(Clipper.ScalePath64(boxPoly, 1000*scale), 1, 2);

        // gcell.addBox(9050, 46500, 1900, 1000, 1, 1);

        // Do we have the cell in the drawing?
        Assert.That(drawing_.findCell("test"), Is.EqualTo(gcell));
       
        g.setDrawing(drawing_);
        g.setValid(true);
        
        if (File.Exists(outDir + "box.gds"))
        {
            File.Delete(outDir + "box.gds");
        }
        gds.gdsWriter gw = new(g, outDir + "box.gds");
        gw.save();
        Assert.That(File.Exists(outDir + "box.gds"), Is.True);

        if (File.Exists(outDir + "box.oas"))
        {
            File.Delete(outDir + "box.oas");
        }
        oasis.oasWriter ow = new(g, outDir + "box.oas");
        ow.save();
        Assert.That(File.Exists(outDir + "box.oas"), Is.True);
    }

    [Test]
    public static void test_cell_export()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(0));

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");
        Assert.That(gcGDS.activeLD, Is.EqualTo(0));

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.That(t1.Count, Is.EqualTo(1));
        Assert.That(t1[0].layer_nr, Is.EqualTo(1));
        Assert.That(t1[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t1[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t1[0].pointarray[0], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t1[0].pointarray[1], Is.EqualTo(new Point64(0,100,0)));
        Assert.That(t1[0].pointarray[2], Is.EqualTo(new Point64(40,100,0)));
        Assert.That(t1[0].pointarray[3], Is.EqualTo(new Point64(40,40,0)));
        Assert.That(t1[0].pointarray[4], Is.EqualTo(new Point64(100,40,0)));
        Assert.That(t1[0].pointarray[5], Is.EqualTo(new Point64(100,0,0)));
        Assert.That(t1[0].pointarray[6], Is.EqualTo(new Point64(0,0,0)));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(1));
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.That(t2.Count, Is.EqualTo(4));
        Assert.That(t2[0].layer_nr, Is.EqualTo(1));
        Assert.That(t2[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[0].pointarray[0], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t2[0].pointarray[1], Is.EqualTo(new Point64(0,100,0)));
        Assert.That(t2[0].pointarray[2], Is.EqualTo(new Point64(40,100,0)));
        Assert.That(t2[0].pointarray[3], Is.EqualTo(new Point64(40,40,0)));
        Assert.That(t2[0].pointarray[4], Is.EqualTo(new Point64(100,40,0)));
        Assert.That(t2[0].pointarray[5], Is.EqualTo(new Point64(100,0,0)));
        Assert.That(t2[0].pointarray[6], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t2[1].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[1].pointarray[0], Is.EqualTo(new Point64(110,0,0)));
        Assert.That(t2[1].pointarray[1], Is.EqualTo(new Point64(110,100,0)));
        Assert.That(t2[1].pointarray[2], Is.EqualTo(new Point64(150,100,0)));
        Assert.That(t2[1].pointarray[3], Is.EqualTo(new Point64(150,40,0)));
        Assert.That(t2[1].pointarray[4], Is.EqualTo(new Point64(210,40,0)));
        Assert.That(t2[1].pointarray[5], Is.EqualTo(new Point64(210,0,0)));
        Assert.That(t2[1].pointarray[6], Is.EqualTo(new Point64(110,0,0)));
        Assert.That(t2[2].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[2].pointarray[0], Is.EqualTo(new Point64(0,110,0)));
        Assert.That(t2[2].pointarray[1], Is.EqualTo(new Point64(0,210,0)));
        Assert.That(t2[2].pointarray[2], Is.EqualTo(new Point64(40,210,0)));
        Assert.That(t2[2].pointarray[3], Is.EqualTo(new Point64(40,150,0)));
        Assert.That(t2[2].pointarray[4], Is.EqualTo(new Point64(100,150,0)));
        Assert.That(t2[2].pointarray[5], Is.EqualTo(new Point64(100,110,0)));
        Assert.That(t2[2].pointarray[6], Is.EqualTo(new Point64(0,110,0)));
        Assert.That(t2[3].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[3].pointarray[0], Is.EqualTo(new Point64(110,110,0)));
        Assert.That(t2[3].pointarray[1], Is.EqualTo(new Point64(110,210,0)));
        Assert.That(t2[3].pointarray[2], Is.EqualTo(new Point64(150,210,0)));
        Assert.That(t2[3].pointarray[3], Is.EqualTo(new Point64(150,150,0)));
        Assert.That(t2[3].pointarray[4], Is.EqualTo(new Point64(210,150,0)));
        Assert.That(t2[3].pointarray[5], Is.EqualTo(new Point64(210,110,0)));
        Assert.That(t2[3].pointarray[6], Is.EqualTo(new Point64(110,110,0)));
        
        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        Assert.That(2, Is.EqualTo(gcGDS.activeStructure));
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        // We converted the drawing, but will use the origin to check the output.
        Assert.That(gcGDS.getDrawing().cellList[gcGDS.activeStructure].elementList[0].isCellrefArray(), Is.True);
        Point64 pos = gcGDS.getDrawing().cellList[gcGDS.activeStructure].elementList[0].getPos();
        Assert.That(t3.Count, Is.EqualTo(16));

        int polIndex = 0;
        int rowPitch = 110;
        int colPitch = 110;
        
        // This nested evaluation has a sequence order of the initial 4 cells in a 2 x 2 placement.
        // This grid is then replicated in its own 2 x 2 placement
        for (int outer_row = 0; outer_row < 2; outer_row++)
        {
            for (int outer_col = 0; outer_col < 2; outer_col++)
            {
                for (int row = 0; row < 2; row++)
                {
                    int y_offset = rowPitch * (row + (outer_row * 2));
                    for (int col = 0; col < 2; col++)
                    {
                        Assert.That(t3[polIndex].layer_nr, Is.EqualTo(1));
                        Assert.That(t3[polIndex].datatype_nr, Is.EqualTo(0));
                        Assert.That(t3[polIndex].pointarray.Count, Is.EqualTo(7));
                        int x_offset = colPitch * (col + (outer_col * 2));
                        Assert.That(t3[polIndex].pointarray[0],
                            Is.EqualTo(new Point64((pos.X + 0) + x_offset, (pos.Y + 0) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[1],
                            Is.EqualTo(new Point64((pos.X + 0) + x_offset, (pos.Y + 100) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[2],
                            Is.EqualTo(new Point64((pos.X + 40) + x_offset, (pos.Y + 100) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[3],
                            Is.EqualTo(new Point64((pos.X + 40) + x_offset, (pos.Y + 40) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[4],
                            Is.EqualTo(new Point64((pos.X + 100) + x_offset, (pos.Y + 40) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[5],
                            Is.EqualTo(new Point64((pos.X + 100) + x_offset, (pos.Y + 0) + y_offset, 0)));
                        Assert.That(t3[polIndex].pointarray[6],
                            Is.EqualTo(new Point64((pos.X + 0) + x_offset, (pos.Y + 0) + y_offset, 0)));
                        polIndex++;
                    }
                }
            }
        }
    }

    [Test]
    public static void consistency_from_oasis()
    {
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(baseDir + "/consistency/triage.oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        string outFile = outDir + "/c3_consistency_from_oas.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcOAS, outFile);
        gw.save();
        Assert.That(File.Exists(outFile), Is.True);

        string outFile2 = outDir + "/c3_consistency_from_oas.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcOAS, outFile2);
        ow.save();
        Assert.That(File.Exists(outFile2), Is.True);
    }

    [Test]
    public static void consistency_from_gds()
    {
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(baseDir + "/consistency/c3_consistency.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        string outFile = outDir + "/c3_consistency_from_gds.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcGDS, outFile);
        gw.save();
        Assert.That(File.Exists(outFile), Is.True);

        string outFile2 = outDir + "/c3_consistency_from_gds.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcGDS, outFile2);
        ow.save();
        Assert.That(File.Exists(outFile2), Is.True);
    }

    [Test]
    public static void test_cell_export_complex()
    {
        // Note that the cell computation is for all layers/dataypes. The list of GCPolygons would need to be filtered separately for the LD of interest.
        string arrayDir = baseDir + "cellrefarray" + Path.DirectorySeparatorChar;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(arrayDir + "L_array_nested_2.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        // Simple cell.
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("r");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(0));

        // Only a single layer datatype.
        gcGDS.activeLD = gcGDS.getActiveStructureLDList().IndexOf("L1D0");
        Assert.That(gcGDS.activeLD, Is.EqualTo(0));

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.That(t1.Count, Is.EqualTo(1));
        Assert.That(t1[0].layer_nr, Is.EqualTo(1));
        Assert.That(t1[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t1[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t1[0].pointarray[0], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t1[0].pointarray[1], Is.EqualTo(new Point64(0,100,0)));
        Assert.That(t1[0].pointarray[2], Is.EqualTo(new Point64(40,100,0)));
        Assert.That(t1[0].pointarray[3], Is.EqualTo(new Point64(40,40,0)));
        Assert.That(t1[0].pointarray[4], Is.EqualTo(new Point64(100,40,0)));
        Assert.That(t1[0].pointarray[5], Is.EqualTo(new Point64(100,0,0)));
        Assert.That(t1[0].pointarray[6], Is.EqualTo(new Point64(0,0,0)));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(1));

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons();
        Assert.That(t2.Count, Is.EqualTo(6));
        Assert.That(t2[0].layer_nr, Is.EqualTo(1));
        Assert.That(t2[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[0].pointarray[0], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t2[0].pointarray[1], Is.EqualTo(new Point64(0,100,0)));
        Assert.That(t2[0].pointarray[2], Is.EqualTo(new Point64(40,100,0)));
        Assert.That(t2[0].pointarray[3], Is.EqualTo(new Point64(40,40,0)));
        Assert.That(t2[0].pointarray[4], Is.EqualTo(new Point64(100,40,0)));
        Assert.That(t2[0].pointarray[5], Is.EqualTo(new Point64(100,0,0)));
        Assert.That(t2[0].pointarray[6], Is.EqualTo(new Point64(0,0,0)));
        Assert.That(t2[1].layer_nr, Is.EqualTo(1));
        Assert.That(t2[1].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[1].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[1].pointarray[0], Is.EqualTo(new Point64(110,0,0)));
        Assert.That(t2[1].pointarray[1], Is.EqualTo(new Point64(110,100,0)));
        Assert.That(t2[1].pointarray[2], Is.EqualTo(new Point64(150,100,0)));
        Assert.That(t2[1].pointarray[3], Is.EqualTo(new Point64(150,40,0)));
        Assert.That(t2[1].pointarray[4], Is.EqualTo(new Point64(210,40,0)));
        Assert.That(t2[1].pointarray[5], Is.EqualTo(new Point64(210,0,0)));
        Assert.That(t2[1].pointarray[6], Is.EqualTo(new Point64(110,0,0)));
        Assert.That(t2[2].layer_nr, Is.EqualTo(1));
        Assert.That(t2[2].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[2].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[2].pointarray[0], Is.EqualTo(new Point64(0,110,0)));
        Assert.That(t2[2].pointarray[1], Is.EqualTo(new Point64(0,210,0)));
        Assert.That(t2[2].pointarray[2], Is.EqualTo(new Point64(40,210,0)));
        Assert.That(t2[2].pointarray[3], Is.EqualTo(new Point64(40,150,0)));
        Assert.That(t2[2].pointarray[4], Is.EqualTo(new Point64(100,150,0)));
        Assert.That(t2[2].pointarray[5], Is.EqualTo(new Point64(100,110,0)));
        Assert.That(t2[2].pointarray[6], Is.EqualTo(new Point64(0,110,0)));
        Assert.That(t2[3].layer_nr, Is.EqualTo(1));
        Assert.That(t2[3].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[3].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[3].pointarray[0], Is.EqualTo(new Point64(110,110,0)));
        Assert.That(t2[3].pointarray[1], Is.EqualTo(new Point64(110,210,0)));
        Assert.That(t2[3].pointarray[2], Is.EqualTo(new Point64(150,210,0)));
        Assert.That(t2[3].pointarray[3], Is.EqualTo(new Point64(150,150,0)));
        Assert.That(t2[3].pointarray[4], Is.EqualTo(new Point64(210,150,0)));
        Assert.That(t2[3].pointarray[5], Is.EqualTo(new Point64(210,110,0)));
        Assert.That(t2[3].pointarray[6], Is.EqualTo(new Point64(110,110,0)));
        Assert.That(t2[4].layer_nr, Is.EqualTo(1));
        Assert.That(t2[4].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[4].pointarray.Count, Is.EqualTo(5));
        Assert.That(t2[4].pointarray[0], Is.EqualTo(new Point64(60,173,0)));
        Assert.That(t2[4].pointarray[1], Is.EqualTo(new Point64(60,201,0)));
        Assert.That(t2[4].pointarray[2], Is.EqualTo(new Point64(83,201,0)));
        Assert.That(t2[4].pointarray[3], Is.EqualTo(new Point64(83,173,0)));
        Assert.That(t2[4].pointarray[4], Is.EqualTo(new Point64(60,173,0)));
        Assert.That(t2[5].layer_nr, Is.EqualTo(2).Within(2));
        Assert.That(t2[5].datatype_nr, Is.EqualTo(0).Within(0));
        Assert.That(t2[5].pointarray.Count, Is.EqualTo(5).Within(5));
        Assert.That(t2[5].pointarray[0], Is.EqualTo(new Point64(52,67,0)));
        Assert.That(t2[5].pointarray[1], Is.EqualTo(new Point64(52,93,0)));
        Assert.That(t2[5].pointarray[2], Is.EqualTo(new Point64(103,93,0)));
        Assert.That(t2[5].pointarray[3], Is.EqualTo(new Point64(103,67,0)));
        Assert.That(t2[5].pointarray[4], Is.EqualTo(new Point64(52,67,0)));

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
        Assert.That(Utils.GetSHA256Hash(t4), Is.EqualTo(Utils.GetSHA256Hash(t3)));
        Assert.That(Utils.GetSHA256Hash(t3a), Is.Not.EqualTo(Utils.GetSHA256Hash(t3)));
        Assert.That(Utils.GetSHA256Hash(t3b), Is.Not.EqualTo(Utils.GetSHA256Hash(t3)));
        Assert.That(Utils.GetSHA256Hash(t4a), Is.EqualTo(Utils.GetSHA256Hash(t3a)));
        Assert.That(Utils.GetSHA256Hash(t4b), Is.EqualTo(Utils.GetSHA256Hash(t3b)));
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
        Assert.That(geo.Count, Is.EqualTo(4));
        
        // We have floats, so this gets a little more awkward.
        Assert.That(Math.Abs(geo[0][0].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][0].y - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][1].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][1].y - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][2].x - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][2].y - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][3].x - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][3].y - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][4].x - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][4].y - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][5].x - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][5].y - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][6].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[0][6].y - 0), Is.LessThanOrEqualTo(1E-13));

        int x_adjust = 110;
        Assert.That(Math.Abs(geo[1][0].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][0].y - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][1].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][1].y - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][2].x - (40 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][2].y - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][3].x - (40 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][3].y - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][4].x - (100 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][4].y - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][5].x - (100 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][5].y - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][6].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[1][6].y - 0), Is.LessThanOrEqualTo(1E-13));

        int y_adjust = 110;
        Assert.That(Math.Abs(geo[2][0].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][0].y - (0 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][1].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][1].y - (100 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][2].x - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][2].y - (100 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][3].x - 40), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][3].y - (40 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][4].x - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][4].y - (40 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][5].x - 100), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][5].y - 0 - y_adjust), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][6].x - 0), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[2][6].y - (0 + y_adjust)), Is.LessThanOrEqualTo(1E-13));

        Assert.That(Math.Abs(geo[3][0].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][1].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][2].x - (40 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][3].x - (40 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][4].x - (100 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][5].x - (100 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][6].x - (0 + x_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][0].y - (0 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][1].y - (100 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][2].y - (100 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][3].y - (40 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][4].y - (40 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][5].y - 0 - y_adjust), Is.LessThanOrEqualTo(1E-13));
        Assert.That(Math.Abs(geo[3][6].y - (0 + y_adjust)), Is.LessThanOrEqualTo(1E-13));
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
        Assert.That(Utils.GetSHA256Hash(geo2), Is.EqualTo("vZoyBBGchu/u+LChif6g6kzH6nGSVj6XHD/GdPwfT7g="));
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

        Assert.That(File.Exists(gdsFile), Is.True);
        
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        // L
        Assert.That(cell_gds.elementList[0].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[0].datatype_nr, Is.EqualTo(0));
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(0));

        // Triangle
        Assert.That(cell_gds.elementList[1].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[1].layer_nr, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[1].datatype_nr, Is.EqualTo(0));
        polys_gds = cell_gds.elementList[1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(0));

        // Pentagram
        Assert.That(cell_gds.elementList[2].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[2].layer_nr, Is.EqualTo(3));
        Assert.That(cell_gds.elementList[2].datatype_nr, Is.EqualTo(0));
        polys_gds = cell_gds.elementList[2].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(6));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(15));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(0));

        // Trapezoid
        Assert.That(cell_gds.elementList[3].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[3].layer_nr, Is.EqualTo(4));
        Assert.That(cell_gds.elementList[3].datatype_nr, Is.EqualTo(0));
        polys_gds = cell_gds.elementList[3].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(15));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(0));

        // Parallelogram
        Assert.That(cell_gds.elementList[4].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[4].layer_nr, Is.EqualTo(5));
        Assert.That(cell_gds.elementList[4].datatype_nr, Is.EqualTo(0));
        polys_gds = cell_gds.elementList[4].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(10));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(0));

        string oasFile = outDir + "simple_polygon.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        // L
        Assert.That(cell_oas.elementList[0].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[0].datatype_nr, Is.EqualTo(0));
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(0));

        // Triangle
        Assert.That(cell_oas.elementList[1].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[1].layer_nr, Is.EqualTo(2));
        Assert.That(cell_oas.elementList[1].datatype_nr, Is.EqualTo(0));
        polys_oas = cell_oas.elementList[1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(0));

        // Pentagram
        Assert.That(cell_oas.elementList[2].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[2].layer_nr, Is.EqualTo(3));
        Assert.That(cell_oas.elementList[2].datatype_nr, Is.EqualTo(0));
        polys_oas = cell_oas.elementList[2].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(6));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(5));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(15));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(5));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(0));

        // Trapezoid
        Assert.That(cell_oas.elementList[3].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[3].layer_nr, Is.EqualTo(4));
        Assert.That(cell_oas.elementList[3].datatype_nr, Is.EqualTo(0));
        polys_oas = cell_oas.elementList[3].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(5));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(15));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(0));

        // Parallelogram
        Assert.That(cell_oas.elementList[4].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[4].layer_nr, Is.EqualTo(5));
        Assert.That(cell_oas.elementList[4].datatype_nr, Is.EqualTo(0));
        polys_oas = cell_oas.elementList[4].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(10));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(0));
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
        // This 5 width becomes an issue with coarse layout grids because it resolves to 3 each side (rather than 2.5), which is observed below.
        // Finer grid resolution will avoid this problem.
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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.That(cell_gds.elementList[0].isPath(), Is.True);
        Assert.That(cell_gds.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[0].datatype_nr, Is.EqualTo(0));
        Path64 path_ = cell_gds.elementList[0].getPath();
        Assert.That(path_.Count, Is.EqualTo(5));
        Assert.That(path_[0].X, Is.EqualTo(0));
        Assert.That(path_[0].Y, Is.EqualTo(0));
        Assert.That(path_[1].X, Is.EqualTo(0));
        Assert.That(path_[1].Y, Is.EqualTo(10));
        Assert.That(path_[2].X, Is.EqualTo(20));
        Assert.That(path_[2].Y, Is.EqualTo(10));
        Assert.That(path_[3].X, Is.EqualTo(20));
        Assert.That(path_[3].Y, Is.EqualTo(40));
        Assert.That(path_[4].X, Is.EqualTo(0));
        Assert.That(path_[4].Y, Is.EqualTo(40));
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(11));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(3));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(3));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(8));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(23));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(8));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(23));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(43));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(43));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(38));
        Assert.That(polys_gds[0].pointarray[7].X, Is.EqualTo(18));
        Assert.That(polys_gds[0].pointarray[7].Y, Is.EqualTo(38));
        Assert.That(polys_gds[0].pointarray[8].X, Is.EqualTo(18));
        Assert.That(polys_gds[0].pointarray[8].Y, Is.EqualTo(13));
        Assert.That(polys_gds[0].pointarray[9].X, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[9].Y, Is.EqualTo(13));
        Assert.That(polys_gds[0].pointarray[10].X, Is.EqualTo(-3));
        Assert.That(polys_gds[0].pointarray[10].Y, Is.EqualTo(-3));

        string oasFile = outDir + "simple_path.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.That(cell_oas.elementList[0].isPath(), Is.True);
        Assert.That(cell_oas.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[0].datatype_nr, Is.EqualTo(0));
        path_ = cell_oas.elementList[0].getPath();
        Assert.That(path_.Count, Is.EqualTo(5));
        Assert.That(path_[0].X, Is.EqualTo(0));
        Assert.That(path_[0].Y, Is.EqualTo(0));
        Assert.That(path_[1].X, Is.EqualTo(0));
        Assert.That(path_[1].Y, Is.EqualTo(10));
        Assert.That(path_[2].X, Is.EqualTo(20));
        Assert.That(path_[2].Y, Is.EqualTo(10));
        Assert.That(path_[3].X, Is.EqualTo(20));
        Assert.That(path_[3].Y, Is.EqualTo(40));
        Assert.That(path_[4].X, Is.EqualTo(0));
        Assert.That(path_[4].Y, Is.EqualTo(40));
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(11));

        // The values here go to 2 rather than 3 (rounded from 2.5) due to OASIS not supporting odd-width paths.
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(2));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(2));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(8));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(22));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(8));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(22));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(42));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(42));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(38));
        Assert.That(polys_oas[0].pointarray[7].X, Is.EqualTo(18));
        Assert.That(polys_oas[0].pointarray[7].Y, Is.EqualTo(38));
        Assert.That(polys_oas[0].pointarray[8].X, Is.EqualTo(18));
        Assert.That(polys_oas[0].pointarray[8].Y, Is.EqualTo(12));
        Assert.That(polys_oas[0].pointarray[9].X, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[9].Y, Is.EqualTo(12));
        Assert.That(polys_oas[0].pointarray[10].X, Is.EqualTo(-2));
        Assert.That(polys_oas[0].pointarray[10].Y, Is.EqualTo(-2));
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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.That(cell_gds.elementList[0].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[0].datatype_nr, Is.EqualTo(0));
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(3600));

        string oasFile = outDir + "simple_circle.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.That(cell_oas.elementList[0].isPolygon(), Is.True);
        Assert.That(cell_oas.elementList[0].layer_nr, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[0].datatype_nr, Is.EqualTo(0));
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(3600));
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

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref1";
        gcell.addCellref();

        gcell.elementList[^1].setPos(new (10, 0));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test"));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref1_mx";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new (10, 0));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test"));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);
        gcell.elementList[^1].setMirrorx();

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

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref2";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new (20, 20));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test2"));
        gcell.elementList[^1].setName("test2");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(2);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref2_mx";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new (20, 20));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test2"));
        gcell.elementList[^1].setName("test2");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(2);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellref.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellref1");
        Assert.That(cell_gds.elementList[^1].isCellref(), Is.True);
        Point64 pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(pos.Y + 0));

        cell_gds = drawing_gds.findCell("test_cellref1_mx");
        Assert.That(cell_gds.elementList[^1].isCellref(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(pos.Y + -40));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(pos.Y + -40));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(pos.Y + -20));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(pos.Y + -20));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(pos.Y + 0));

        cell_gds = drawing_gds.findCell("test_cellref2");
        Assert.That(cell_gds.elementList[^1].isCellref(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(20));
        Assert.That(pos.Y, Is.EqualTo(20));
        Assert.That(cell_gds.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(pos.X - (2 * 30)));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(pos.X - (2 * 30)));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(pos.Y));

        cell_gds = drawing_gds.findCell("test_cellref2_mx");
        Assert.That(cell_gds.elementList[^1].isCellref(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(20));
        Assert.That(pos.Y, Is.EqualTo(20));
        Assert.That(cell_gds.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(270));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(pos.X - 60));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(pos.X - 60));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(pos.Y - 20));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(pos.Y - 20));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(pos.Y - 40));
        Assert.That(polys_gds[0].pointarray[5].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[5].Y, Is.EqualTo(pos.Y - 40));
        Assert.That(polys_gds[0].pointarray[6].X, Is.EqualTo(pos.X));
        Assert.That(polys_gds[0].pointarray[6].Y, Is.EqualTo(pos.Y));

        string oasFile = outDir + "simple_cellref.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellref1");
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(pos.Y + 0));

        cell_oas = drawing_oas.findCell("test_cellref1_mx");
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(pos.Y + -40));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(pos.Y + -40));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(pos.X + 20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(pos.Y + -20));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(pos.Y + -20));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(pos.X + 40));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(pos.Y + 0));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(pos.X + 0));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(pos.Y + 0));
        
        cell_oas = drawing_gds.findCell("test_cellref2");
        Assert.That(cell_oas.elementList[^1].isCellref(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(20));
        Assert.That(pos.Y, Is.EqualTo(20));
        Assert.That(cell_oas.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(pos.X - (2 * 30)));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(pos.X - (2 * 30)));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(pos.Y + 20));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(pos.Y + 40));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(pos.Y));
        
        cell_oas = drawing_gds.findCell("test_cellref2_mx");
        Assert.That(cell_oas.elementList[^1].isCellref(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(20));
        Assert.That(pos.Y, Is.EqualTo(20));
        Assert.That(cell_oas.elementList[^1].getScale(), Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(270));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(polys_oas[0].pointarray[0].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[0].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_oas[0].pointarray[1].X, Is.EqualTo(pos.X - 60));
        Assert.That(polys_oas[0].pointarray[1].Y, Is.EqualTo(pos.Y));
        Assert.That(polys_oas[0].pointarray[2].X, Is.EqualTo(pos.X - 60));
        Assert.That(polys_oas[0].pointarray[2].Y, Is.EqualTo(pos.Y - 20));
        Assert.That(polys_oas[0].pointarray[3].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_oas[0].pointarray[3].Y, Is.EqualTo(pos.Y - 20));
        Assert.That(polys_oas[0].pointarray[4].X, Is.EqualTo(pos.X - 20));
        Assert.That(polys_oas[0].pointarray[4].Y, Is.EqualTo(pos.Y - 40));
        Assert.That(polys_oas[0].pointarray[5].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[5].Y, Is.EqualTo(pos.Y - 40));
        Assert.That(polys_oas[0].pointarray[6].X, Is.EqualTo(pos.X));
        Assert.That(polys_oas[0].pointarray[6].Y, Is.EqualTo(pos.Y));
    }

    [Test]
    public static void defineAndWrite_CellrefArray1()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray1";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (0, 0));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray1.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray1");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(0));
        Assert.That(pos.Y, Is.EqualTo(0));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray1.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray1");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(0));
        Assert.That(pos.Y, Is.EqualTo(0));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray2()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray2";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);
        
        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray2.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray2");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray2.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray2");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray3()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray3";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(1);

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray3.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray3");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(30 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
        
        string oasFile = outDir + "simple_cellrefarray3.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray3");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));

        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(30 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(0 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(0 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray4()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray4";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);
        
        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray4.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray4");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        Point64 pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        Point64 count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        Point64 col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        Point64 row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        double scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        int polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray4.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray4");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));

        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray5()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray5";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(2);

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray5.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray5");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        Point64 pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        Point64 count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        Point64 col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        Point64 row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        double scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        int polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(-30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(-30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray5.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray5");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));

        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(-30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(-30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(-10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray1mx()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray1_mx";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (0, 0));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray1_mx.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray1_mx");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(0));
        Assert.That(pos.Y, Is.EqualTo(0));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray1_mx.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray1_mx");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(0));
        Assert.That(pos.Y, Is.EqualTo(0));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray2mx()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray2_mx";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray2_mx.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;
        
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray2_mx");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray2_mx.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray2_mx");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * col_pitch.X)))));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * row_pitch.Y)))));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray3mx()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray3_mx";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(1);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray3_mx.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray3_mx");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(30 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(20 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(20 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray3_mx.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray3_mx");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));

        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(30 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(20 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(20 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray4mx()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray4_mx";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray4_mx.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray4_mx");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray4_mx.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray4_mx");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void defineAndWrite_CellrefArray5mx()
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
        array[1] = new (0, 80);
        array[2] = new (100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray5_mx";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new (10, 20));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(2);
        gcell.elementList[^1].setMirrorx();

        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "simple_cellrefarray5_mx.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        Point64 pos, row_pitch, col_pitch, count;
        double scale;
        int polyIndex = 0;

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test_cellrefarray5_mx");
        Assert.That(cell_gds.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_gds.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_gds.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(50 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(50 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(40 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }

        string oasFile = outDir + "simple_cellrefarray5_mx.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test_cellrefarray5_mx");
        Assert.That(cell_oas.elementList[^1].isCellrefArray(), Is.True);
        pos = cell_oas.elementList[^1].getPos();
        Assert.That(pos.X, Is.EqualTo(10));
        Assert.That(pos.Y, Is.EqualTo(20));
        count = cell_oas.elementList[^1].getCount();
        Assert.That(count.X, Is.EqualTo(4));
        Assert.That(count.Y, Is.EqualTo(4));
        col_pitch = cell_oas.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(0 / count.Y));
        row_pitch = cell_oas.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_oas.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_oas.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_oas.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(16));
        
        polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_oas[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_oas[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[1].X, Is.EqualTo(50 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[2].X, Is.EqualTo(50 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(40 +  (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[3].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[4].X, Is.EqualTo(30 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * row_pitch.Y)));
                Assert.That(polys_oas[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * col_pitch.X)));
                Assert.That(polys_oas[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * row_pitch.Y)));
                polyIndex++;
            }
        }
    }

    [Test]
    public static void read_cblock_oasis()
    {
        // Bring in our reference file for comparison sake. This gets complicated.
        string oasFile_ref = baseDir + "compression_test.oas";
        GeoCoreHandler gH_OAS_ref = new();
        gH_OAS_ref.updateGeoCoreHandler(oasFile_ref, GeoCore.fileType.oasis);
        GeoCore gcOAS_ref = gH_OAS_ref.getGeo();
        Assert.That(gcOAS_ref.isValid(), Is.True);
        /*
        GCDrawingfield drawing_oas_ref = gcOAS_ref.getDrawing();
        GCCell cell_oas_ref = drawing_oas_ref.findCell("test_cellrefarray1");
        List<GCPolygon> polys_oas_ref = cell_oas_ref.elementList[^1].convertToPolygons();
        */
    }

    // FIXME : This needs the test conditions added to ensure we can read what we write, etc.
    [Test]
    public static void defineAndWrite_CellrefArray_irregular()
    {
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

        Path64 array = Helper.initedPath64(3);
        array[0] = new (0, 0);
        array[1] = new (30, 50);
        array[2] = new (90, 70);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray_irregular";
        gcell.addCellref(drawing_.findCell("test"), new (0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 2, 2);
        gcell.elementList[^1].setPos(new (0, 0));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(1);
        
        g.setDrawing(drawing_);
        g.setValid(true);

        string gdsFile = outDir + "cellref_array_irregular.gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(g, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        string oasFile = outDir + "cellref_array_irregular.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
    }
    
    // Need tests for the interaction of array position and so on as well, e.g. the initial array 0 entry, pos, and so on.

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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.That(cell_gds.elementList[^1].isText(), Is.True);

        string oasFile = outDir + "simple_text.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.That(cell_oas.elementList[^1].isText(), Is.True);
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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.That(cell_gds.elementList[^1].isPolygon(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(5));

        string oasFile = outDir + "simple_ctrapezoid.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.That( cell_oas.elementList[^1].isPolygon(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(5));
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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("test");

        Assert.That(cell_gds.elementList[^1].isPolygon(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
        Assert.That(polys_gds.Count, Is.EqualTo(1));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(5));

        string oasFile = outDir + "simple_trapezoid.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("test");

        Assert.That(cell_oas.elementList[^1].isPolygon(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
        Assert.That(polys_oas.Count, Is.EqualTo(1));
        Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(5));
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
        Assert.That(File.Exists(gdsFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);
        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        for (int i = 0; i < edge * edge; i++)
        {
            GCCell cell_gds = drawing_gds.findCell("test" + i);

            Assert.That(cell_gds.elementList[^1].isPolygon(), Is.True);
            List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons();
            Assert.That(polys_gds.Count, Is.EqualTo(1));
            Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(7));
        }

        string oasFile = outDir + edge + "_cellref.oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(g, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);
        
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);
        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        for (int i = 0; i < edge * edge; i++)
        {
            GCCell cell_oas = drawing_oas.findCell("test" + i);

            Assert.That(cell_oas.elementList[^1].isPolygon(), Is.True);
            List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons();
            Assert.That(polys_oas.Count, Is.EqualTo(1));
            Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        }
    }

    [Test]
    public static void move_box_test()
    {
        int dimension = 10;
        string filename = "move_box_" + 10 + "_" + 20;
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
        Assert.That(2018, Is.EqualTo(drawing_.accyear));

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

        gcell.addBox(0 * scale, 0 * scale, dimension * scale, dimension * scale, 1, 1);
        GCElement content = gcell.elementList[^1];

        Assert.That(content.isBox(), Is.True);
        Assert.That(content.getPos().X, Is.EqualTo(0));
        Assert.That(content.getPos().Y, Is.EqualTo(0));
        List<GCPolygon> polys = content.convertToPolygons();
        Int64 minx = Int64.MaxValue;
        Int64 miny = Int64.MaxValue;
        foreach (GCPolygon pl in polys)
        {
            Int64 tmp = pl.pointarray.MinBy(p => p.X).X;
            if (tmp < minx)
            {
                minx = tmp;
            }
            tmp = pl.pointarray.MinBy(p => p.Y).Y;
            if (tmp < miny)
            {
                miny = tmp;
            }
        }
        Assert.That(minx, Is.EqualTo(0));
        Assert.That(miny, Is.EqualTo(0));

        int move_x = 10;
        int move_y = 20;
        content.move(new (move_x * scale, move_y * scale));
        Assert.That(content.getPos().X, Is.EqualTo(move_x * scale));
        Assert.That(content.getPos().Y, Is.EqualTo(move_y * scale));
        minx = Int64.MaxValue;
        miny = Int64.MaxValue;
        polys = content.convertToPolygons();
        foreach (GCPolygon pl in polys)
        {
            Int64 tmp = pl.pointarray.MinBy(p => p.X).X;
            if (tmp < minx)
            {
                minx = tmp;
            }
            tmp = pl.pointarray.MinBy(p => p.Y).Y;
            if (tmp < miny)
            {
                miny = tmp;
            }
        }
        Assert.That(minx, Is.EqualTo(move_x * scale));
        Assert.That(miny, Is.EqualTo(move_y * scale));
        g.setDrawing(drawing_);
        g.setValid(true);

        if (File.Exists(outDir + "/" + filename + ".gds"))
        {
            File.Delete(outDir + "/" + filename + ".gds");
        }

        gds.gdsWriter gw = new(g, outDir + "/" + filename + ".gds");
        gw.save();
        Assert.That(File.Exists(outDir + "/" + filename + ".gds"), Is.True);

        if (File.Exists(outDir + "/" + filename + ".oas"))
        {
            File.Delete(outDir + "/" + filename + ".oas");
        }

        oasis.oasWriter ow = new(g, outDir + "/" + filename + ".oas");
        ow.save();
        Assert.That(File.Exists(outDir + "/" + filename + ".oas"), Is.True);

        // Load the files in to see what we have.

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(outDir + "/" + filename + ".gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        drawing_gds.databaseunits = 1000 * scale;
        drawing_gds.userunits = 0.001 / scale;
        GCCell cell_gds = drawing_gds.findCell("test");
        int elementCount = cell_gds.elementList.Count;
        Assert.That(elementCount, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].isBox(), Is.True);
        Assert.That(cell_gds.elementList[^1].getPos().X, Is.EqualTo(1000));
        Assert.That(cell_gds.elementList[^1].getPos().Y, Is.EqualTo(2000));

        string out_filename = outDir + "/" + filename + "_resave_from_gds.gds";
        save_gdsii(gcGDS, out_filename);

        out_filename = outDir + "/" + filename + "_resave_from_gds.oas";
        save_oasis(gcGDS, out_filename);

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(outDir + "/" + filename + ".oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        drawing_oas.databaseunits = 1000 * scale;
        drawing_oas.userunits = 0.001 / scale;
        GCCell cell_oas = drawing_oas.findCell("test");
        elementCount = cell_oas.elementList.Count;
        Assert.That(elementCount, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].isBox(), Is.True);
        Assert.That(cell_oas.elementList[^1].getPos().X, Is.EqualTo(1000));
        Assert.That(cell_oas.elementList[^1].getPos().Y, Is.EqualTo(2000));

        out_filename = outDir + "/" + filename + "_resave_from_oas.gds";
        save_gdsii(gcOAS, out_filename);

        out_filename = outDir + "/" + filename + "_resave_from_oas.oas";
        save_oasis(gcOAS, out_filename);
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
            Assert.That(File.Exists(gdsFile), Is.True);

            string oasFile = outDir + "L" + layer + "D" + datatype + "_box.oas";
            if (File.Exists(oasFile))
            {
                File.Delete(oasFile);
            }
            oasis.oasWriter ow = new(g, oasFile);
            ow.save();
            Assert.That(File.Exists(oasFile), Is.True);
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
            Assert.That(File.Exists(gdsFile), Is.True);

            string oasFile = outDir + "L" + layer + "D" + datatype + "_box.oas";
            if (File.Exists(oasFile))
            {
                File.Delete(oasFile);
            }
            oasis.oasWriter ow = new(g, oasFile);
            ow.save();
            Assert.That(File.Exists(oasFile), Is.True);
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
        Assert.That(File.Exists(out_filename), Is.True);
    }

    private static void save_oasis(GeoCore gc, string out_filename)
    {
        if (File.Exists(out_filename))
        {
            File.Delete(out_filename);
        }
        oasis.oasWriter ow = new(gc, out_filename);
        ow.save();
        Assert.That(File.Exists(out_filename), Is.True);
    }
}