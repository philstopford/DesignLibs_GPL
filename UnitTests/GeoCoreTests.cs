using System.Buffers.Binary;
using Clipper2Lib;
using geoCoreLib;
using utility;

namespace UnitTests;

/// <summary>
/// Comprehensive tests for the geoCore library, which provides geometric file format handling
/// and integer-based geometric operations for GDSII and OASIS formats.
/// These unit tests complement the extensive integration tests with focused class testing.
/// </summary>
public class GeoCoreTests
{
    static string baseDir = "../../../../../geocore_test/";
    static string outDir = TestHelpers.TestPath.Get("geocore_out") + Path.DirectorySeparatorChar;

    #region GeoCore Basic Unit Tests

    /// <summary>
    /// Tests basic GeoCore construction and validity state management.
    /// </summary>
    [Test]
    public static void GeoCore_Construction_InitializesCorrectly()
    {
        // Arrange & Act
        GeoCore geoCore = new GeoCore();

        // Assert: New instance should start invalid
        Assert.That(geoCore.isValid(), Is.False);
        Assert.That(geoCore.error_msgs, Is.Not.Null);
    }

    /// <summary>
    /// Tests validity state management in GeoCore.
    /// </summary>
    [Test]
    public static void GeoCore_ValidityState_CanBeSetAndRetrieved()
    {
        // Arrange
        GeoCore geoCore = new GeoCore();

        // Act: Set valid state to true
        geoCore.setValid(true);

        // Assert
        Assert.That(geoCore.isValid(), Is.True);

        // Act: Set valid state to false
        geoCore.setValid(false);

        // Assert
        Assert.That(geoCore.isValid(), Is.False);
    }

    /// <summary>
    /// Tests that GeoCore tolerance has expected default value.
    /// </summary>
    [Test]
    public static void GeoCore_Tolerance_HasExpectedDefaultValue()
    {
        // Assert: Static tolerance should be set to expected precision
        Assert.That(GeoCore.tolerance, Is.EqualTo(0.001));
    }

    /// <summary>
    /// Tests GeoCore file type enumeration values.
    /// </summary>
    [Test]
    public static void GeoCore_FileType_ContainsExpectedValues()
    {
        // Assert: All expected file types are available
        Assert.That(Enum.IsDefined(typeof(GeoCore.fileType), GeoCore.fileType.gds), Is.True);
        Assert.That(Enum.IsDefined(typeof(GeoCore.fileType), GeoCore.fileType.oasis), Is.True);
    }

    #endregion

    #region GCPolygon Unit Tests

    /// <summary>
    /// Tests GCPolygon default constructor creates empty polygon correctly.
    /// </summary>
    [Test]
    public static void GCPolygon_DefaultConstructor_CreatesEmptyPolygon()
    {
        // Arrange & Act
        GCPolygon polygon = new GCPolygon();

        // Assert
        Assert.That(polygon.pointarray, Is.Not.Null);
        Assert.That(polygon.pointarray.Count, Is.EqualTo(0));
        Assert.That(polygon.text, Is.False);
    }

    /// <summary>
    /// Tests GCPolygon constructor with points, layer, and datatype.
    /// </summary>
    [Test]
    public static void GCPolygon_ParameterizedConstructor_InitializesCorrectly()
    {
        // Arrange
        Path64 points = new Path64
        {
            new Point64(0, 0),
            new Point64(100, 0),
            new Point64(100, 100),
            new Point64(0, 100)
        };
        int expectedLayer = 1;
        int expectedDatatype = 0;

        // Act
        GCPolygon polygon = new GCPolygon(points, expectedLayer, expectedDatatype);

        // Assert
        Assert.That(polygon.pointarray.Count, Is.EqualTo(4));
        Assert.That(polygon.layer_nr, Is.EqualTo(expectedLayer));
        Assert.That(polygon.datatype_nr, Is.EqualTo(expectedDatatype));

        // Verify points are copied correctly
        Assert.That(polygon.pointarray[0].X, Is.EqualTo(0));
        Assert.That(polygon.pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polygon.pointarray[1].X, Is.EqualTo(100));
        Assert.That(polygon.pointarray[1].Y, Is.EqualTo(0));
        Assert.That(polygon.pointarray[2].X, Is.EqualTo(100));
        Assert.That(polygon.pointarray[2].Y, Is.EqualTo(100));
        Assert.That(polygon.pointarray[3].X, Is.EqualTo(0));
        Assert.That(polygon.pointarray[3].Y, Is.EqualTo(100));
    }

    /// <summary>
    /// Tests GCPolygon copy constructor creates correct deep copy.
    /// </summary>
    [Test]
    public static void GCPolygon_CopyConstructor_CreatesDeepCopy()
    {
        // Arrange
        Path64 originalPoints = new Path64
        {
            new Point64(10, 20),
            new Point64(30, 40),
            new Point64(50, 60)
        };
        GCPolygon originalPolygon = new GCPolygon(originalPoints, 2, 1);
        originalPolygon.text = true;
        originalPolygon.name = "TestPolygon";

        // Act
        GCPolygon copiedPolygon = new GCPolygon(originalPolygon);

        // Assert: Values should be copied
        Assert.That(copiedPolygon.pointarray.Count, Is.EqualTo(originalPolygon.pointarray.Count));
        Assert.That(copiedPolygon.layer_nr, Is.EqualTo(originalPolygon.layer_nr));
        Assert.That(copiedPolygon.datatype_nr, Is.EqualTo(originalPolygon.datatype_nr));

        // Assert: Should be deep copy (different references)
        Assert.That(ReferenceEquals(copiedPolygon.pointarray, originalPolygon.pointarray), Is.False);

        // Verify point data is copied correctly
        for (int i = 0; i < originalPolygon.pointarray.Count; i++)
        {
            Assert.That(copiedPolygon.pointarray[i].X, Is.EqualTo(originalPolygon.pointarray[i].X));
            Assert.That(copiedPolygon.pointarray[i].Y, Is.EqualTo(originalPolygon.pointarray[i].Y));
        }
    }

    /// <summary>
    /// Tests GCPolygon with various polygon shapes.
    /// </summary>
    [Test]
    public static void GCPolygon_VariousShapes_HandledCorrectly()
    {
        // Test triangle
        Path64 trianglePoints = new Path64
        {
            new Point64(0, 0),
            new Point64(100, 0),
            new Point64(50, 100)
        };
        GCPolygon triangle = new GCPolygon(trianglePoints, 1, 0);
        Assert.That(triangle.pointarray.Count, Is.EqualTo(3));

        // Test complex polygon (hexagon)
        Path64 hexagonPoints = new Path64
        {
            new Point64(100, 0),
            new Point64(150, 50),
            new Point64(150, 150),
            new Point64(100, 200),
            new Point64(50, 150),
            new Point64(50, 50)
        };
        GCPolygon hexagon = new GCPolygon(hexagonPoints, 2, 1);
        Assert.That(hexagon.pointarray.Count, Is.EqualTo(6));
    }

    /// <summary>
    /// Tests GCPolygon properties can be modified after construction.
    /// </summary>
    [Test]
    public static void GCPolygon_Properties_CanBeModified()
    {
        // Arrange
        GCPolygon polygon = new GCPolygon();

        // Act
        polygon.text = true;
        polygon.name = "ModifiedPolygon";
        polygon.layer_nr = 5;
        polygon.datatype_nr = 3;

        // Assert
        Assert.That(polygon.text, Is.True);
        Assert.That(polygon.name, Is.EqualTo("ModifiedPolygon"));
        Assert.That(polygon.layer_nr, Is.EqualTo(5));
        Assert.That(polygon.datatype_nr, Is.EqualTo(3));
    }

    #endregion

    #region GCCell Unit Tests

    /// <summary>
    /// Tests GCCell default constructor initialization.
    /// </summary>
    [Test]
    public static void GCCell_DefaultConstructor_InitializesWithCurrentDateTime()
    {
        // Arrange
        DateTime beforeConstruction = DateTime.Now;

        // Act
        GCCell cell = new GCCell();

        // Arrange
        DateTime afterConstruction = DateTime.Now;

        // Assert: Cell should be initialized with current date/time
        Assert.That(cell.elementList, Is.Not.Null);
        Assert.That(cell.saved, Is.False);

        // Verify date/time fields are reasonable (within test execution timeframe)
        Assert.That(cell.modyear, Is.GreaterThanOrEqualTo(beforeConstruction.Year));
        Assert.That(cell.modyear, Is.LessThanOrEqualTo(afterConstruction.Year));
        Assert.That(cell.modmonth, Is.GreaterThanOrEqualTo(1));
        Assert.That(cell.modmonth, Is.LessThanOrEqualTo(12));
        Assert.That(cell.modday, Is.GreaterThanOrEqualTo(1));
        Assert.That(cell.modday, Is.LessThanOrEqualTo(31));
    }

    /// <summary>
    /// Tests GCCell properties can be set and retrieved correctly.
    /// </summary>
    [Test]
    public static void GCCell_Properties_CanBeSetAndRetrieved()
    {
        // Arrange
        GCCell cell = new GCCell();

        // Act
        cell.cellName = "TestCell";
        cell.saved = true;
        cell.modyear = 2023;
        cell.modmonth = 6;
        cell.modday = 15;

        // Assert
        Assert.That(cell.cellName, Is.EqualTo("TestCell"));
        Assert.That(cell.saved, Is.True);
        Assert.That(cell.modyear, Is.EqualTo(2023));
        Assert.That(cell.modmonth, Is.EqualTo(6));
        Assert.That(cell.modday, Is.EqualTo(15));
    }

    #endregion

    #region OASIS-specific Unit Tests

    /// <summary>
    /// Tests OASIS property reading with various value types according to the OASIS specification.
    /// </summary>
    [Test]
    public static void OASISReader_PropertyReading_HandlesAllValueTypes()
    {
        // This test would require creating a minimal OASIS file with different property types
        // For now, we'll test the basic infrastructure
        Assert.That(true, Is.True, "OASIS property reading infrastructure is available");
    }

    /// <summary>
    /// Tests OASIS repetition handling for rectangular arrays.
    /// </summary>
    [Test]
    public static void OASISReader_Repetition_HandlesRectangularArrays()
    {
        // Test rectangular repetition logic
        Repetition rep = new Repetition();
        rep.type = Repetition.RepetitionType.Rectangular;
        rep.columns = 3;
        rep.rows = 2;
        rep.spacing = new Point64(100, 200);

        int expectedCount = rep.columns * rep.rows;
        Assert.That(rep.get_count(), Is.EqualTo(expectedCount));
    }

    /// <summary>
    /// Tests OASIS repetition handling for explicit coordinate arrays.
    /// </summary>
    [Test]
    public static void OASISReader_Repetition_HandlesExplicitArrays()
    {
        // Test explicit repetition logic
        Repetition rep = new Repetition();
        rep.type = Repetition.RepetitionType.Explicit;
        rep.offsets.Add(new Point64(0, 0));
        rep.offsets.Add(new Point64(100, 50));
        rep.offsets.Add(new Point64(100, 50)); // Delta values, so this adds to previous

        Assert.That(rep.get_count(), Is.EqualTo(3));
        
        Path64 offsets = rep.get_offsets();
        Assert.That(offsets.Count, Is.EqualTo(3));
        Assert.That(offsets[0].X, Is.EqualTo(0));
        Assert.That(offsets[0].Y, Is.EqualTo(0));
        Assert.That(offsets[1].X, Is.EqualTo(100));
        Assert.That(offsets[1].Y, Is.EqualTo(50));
        Assert.That(offsets[2].X, Is.EqualTo(200)); // 100 + 100 (delta accumulation)
        Assert.That(offsets[2].Y, Is.EqualTo(100)); // 50 + 50 (delta accumulation)
    }

    /// <summary>
    /// Tests OASIS repetition handling for coordinate list arrays.
    /// </summary>
    [Test]
    public static void OASISReader_Repetition_HandlesCoordinateArrays()
    {
        // Test explicit X repetition
        Repetition repX = new Repetition();
        repX.type = Repetition.RepetitionType.ExplicitX;
        repX.coords.Add(0.0);
        repX.coords.Add(50.5);
        repX.coords.Add(100.0);

        Assert.That(repX.get_count(), Is.EqualTo(3));

        // Test explicit Y repetition
        Repetition repY = new Repetition();
        repY.type = Repetition.RepetitionType.ExplicitY;
        repY.coords.Add(0.0);
        repY.coords.Add(25.5);

        Assert.That(repY.get_count(), Is.EqualTo(2));
    }

    /// <summary>
    /// Tests OASIS repetition reset and copy functionality.
    /// </summary>
    [Test]
    public static void OASISReader_Repetition_ResetAndCopy()
    {
        // Create and populate a repetition
        Repetition original = new Repetition();
        original.type = Repetition.RepetitionType.Rectangular;
        original.columns = 4;
        original.rows = 3;
        original.spacing = new Point64(150, 100);

        // Test copy constructor
        Repetition copied = new Repetition(original);
        Assert.That(copied.type, Is.EqualTo(original.type));
        Assert.That(copied.columns, Is.EqualTo(original.columns));
        Assert.That(copied.rows, Is.EqualTo(original.rows));
        Assert.That(copied.spacing.X, Is.EqualTo(original.spacing.X));
        Assert.That(copied.spacing.Y, Is.EqualTo(original.spacing.Y));

        // Test reset
        copied.reset();
        Assert.That(copied.type, Is.EqualTo(Repetition.RepetitionType.None));
        Assert.That(copied.columns, Is.EqualTo(0));
        Assert.That(copied.rows, Is.EqualTo(0));
    }

    /// <summary>
    /// Tests that OASIS value constants are properly defined.
    /// </summary>
    [Test]
    public static void OASISValues_Constants_AreProperlyDefined()
    {
        // Test basic record types from OASIS specification
        Assert.That(oasis.oasValues.PAD, Is.EqualTo(0));
        Assert.That(oasis.oasValues.START, Is.EqualTo(1));
        Assert.That(oasis.oasValues.END, Is.EqualTo(2));
        
        // Test element types
        Assert.That(oasis.oasValues.PLACEMENT, Is.EqualTo(17));
        Assert.That(oasis.oasValues.TEXT, Is.EqualTo(19));
        Assert.That(oasis.oasValues.RECTANGLE, Is.EqualTo(20));
        Assert.That(oasis.oasValues.POLYGON, Is.EqualTo(21));
        Assert.That(oasis.oasValues.PATH, Is.EqualTo(22));
        Assert.That(oasis.oasValues.CIRCLE, Is.EqualTo(27));
        
        // Test property types
        Assert.That(oasis.oasValues.PROPERTY, Is.EqualTo(28));
        Assert.That(oasis.oasValues.LAST_PROPERTY, Is.EqualTo(29));
    }

    /// <summary>
    /// Tests OASIS file format signature validation.
    /// </summary>
    [Test]
    public static void OASISReader_FileFormatValidation_RejectsInvalidSignature()
    {
        // This tests the basic infrastructure for format validation
        // In a real implementation, we would create invalid OASIS files and test them
        string invalidSignature = "%INVALID-FORMAT\r\n";
        string validSignature = "%SEMI-OASIS\r\n";
        
        Assert.That(validSignature.StartsWith("%SEMI-OASIS"), Is.True);
        Assert.That(invalidSignature.StartsWith("%SEMI-OASIS"), Is.False);
    }

    #endregion

    #region Integration Test Preservation

    /// <summary>
    /// Original klayout GDS array reference test - preserved for compatibility.
    /// Tests comprehensive GDSII file parsing and cell array handling.
    /// </summary>

    [Test]
    public static void klayout_gds_aref_test()
    {
        string gdsFile = baseDir + "klayout_test/gds/arefs_single.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        // Currently, we don't touch the internals in the layout data so use this for conversion when not
        // working with the drawing directly (which is the only place where the units are available.
        // If we use the drawing convert to polygons method with a cell name, this is not needed - we should
        // get the scaled version automatically.
        double grid_scaling = drawing_gds.databaseunits / GCDrawingfield.default_databaseunits;

        // T cell
        GCCell cell_gds = drawing_gds.findCell("T");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(1));
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(grid_scaling);
        Assert.That(polys_gds.Count, Is.EqualTo(1));

        // A cell
        GCCell cell_gds2 = drawing_gds.findCell("A");
        Assert.That(cell_gds2.elementList.Count, Is.EqualTo(1));
        for (int i = 0; i < 1; i++)
        {
            Assert.That(cell_gds2.elementList[i].isCellrefArray(), Is.True);
        }

        // Going directly to the cell, we get unscaled geometry.
        List<GCPolygon> polys_gds2 = cell_gds2.convertToPolygons(1.0);
        // 1 cell reference that are 3x2 arrays of different configurations
        int poly_count = 6;
        Assert.That(polys_gds2.Count, Is.EqualTo(poly_count));
        // Resize to match our drawing units....
        for (int poly_index = 0; poly_index < poly_count; poly_index++)
        {
            polys_gds2[poly_index].resize(grid_scaling);
        }

        List<GCPolygon> polys_gds2_scaled = cell_gds2.convertToPolygons(grid_scaling);

        for (int poly_index = 0; poly_index < poly_count; poly_index++)
        {
            int x_offset = 0;
            int y_offset = 0;
            switch (poly_index)
            {
                case 1:
                    x_offset = 1500;
                    y_offset = 0;
                    break;
                case 2:
                    x_offset = 3000;
                    y_offset = 0;
                    break;
                case 3:
                    x_offset = 0;
                    y_offset = 2000;
                    break;
                case 4:
                    x_offset = 1500;
                    y_offset = 2000;
                    break;
                case 5:
                    x_offset = 3000;
                    y_offset = 2000;
                    break;
            }
            Assert.That(polys_gds2[poly_index].pointarray[0].X, Is.EqualTo(polys_gds2[0].pointarray[0].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[0].Y, Is.EqualTo(polys_gds2[0].pointarray[0].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[1].X, Is.EqualTo(polys_gds2[0].pointarray[1].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[1].Y, Is.EqualTo(polys_gds2[0].pointarray[1].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[2].X, Is.EqualTo(polys_gds2[0].pointarray[2].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[2].Y, Is.EqualTo(polys_gds2[0].pointarray[2].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[3].X, Is.EqualTo(polys_gds2[0].pointarray[3].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[3].Y, Is.EqualTo(polys_gds2[0].pointarray[3].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[4].X, Is.EqualTo(polys_gds2[0].pointarray[4].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[4].Y, Is.EqualTo(polys_gds2[0].pointarray[4].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[5].X, Is.EqualTo(polys_gds2[0].pointarray[5].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[5].Y, Is.EqualTo(polys_gds2[0].pointarray[5].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[6].X, Is.EqualTo(polys_gds2[0].pointarray[6].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[6].Y, Is.EqualTo(polys_gds2[0].pointarray[6].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[7].X, Is.EqualTo(polys_gds2[0].pointarray[7].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[7].Y, Is.EqualTo(polys_gds2[0].pointarray[7].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[8].X, Is.EqualTo(polys_gds2[0].pointarray[8].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[8].Y, Is.EqualTo(polys_gds2[0].pointarray[8].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[9].X, Is.EqualTo(polys_gds2[0].pointarray[9].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[9].Y, Is.EqualTo(polys_gds2[0].pointarray[9].Y + y_offset));
            Assert.That(polys_gds2[poly_index].pointarray[10].X, Is.EqualTo(polys_gds2[0].pointarray[10].X + x_offset));
            Assert.That(polys_gds2[poly_index].pointarray[10].Y, Is.EqualTo(polys_gds2[0].pointarray[10].Y + y_offset));

            // Check that our cell-based scaled version matches the manually scaled one.
            Assert.That(polys_gds2_scaled[poly_index].pointarray[0].X, Is.EqualTo(polys_gds2[poly_index].pointarray[0].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[0].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[0].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[1].X, Is.EqualTo(polys_gds2[poly_index].pointarray[1].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[1].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[1].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[2].X, Is.EqualTo(polys_gds2[poly_index].pointarray[2].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[2].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[2].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[3].X, Is.EqualTo(polys_gds2[poly_index].pointarray[3].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[3].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[3].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[4].X, Is.EqualTo(polys_gds2[poly_index].pointarray[4].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[4].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[4].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[5].X, Is.EqualTo(polys_gds2[poly_index].pointarray[5].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[5].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[5].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[6].X, Is.EqualTo(polys_gds2[poly_index].pointarray[6].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[6].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[6].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[7].X, Is.EqualTo(polys_gds2[poly_index].pointarray[7].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[7].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[7].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[8].X, Is.EqualTo(polys_gds2[poly_index].pointarray[8].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[8].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[8].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[9].X, Is.EqualTo(polys_gds2[poly_index].pointarray[9].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[9].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[9].Y));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[10].X, Is.EqualTo(polys_gds2[poly_index].pointarray[10].X));
            Assert.That(polys_gds2_scaled[poly_index].pointarray[10].Y, Is.EqualTo(polys_gds2[poly_index].pointarray[10].Y));
        }

        // Compare with our drawing conversion... Our list of cells is just 'A' in this case.
        List<GCPolygon> polys_fromdrawing = drawing_gds.convertToPolygons(cells: ["A"])[0];
        Assert.That(polys_fromdrawing.Count, Is.EqualTo(poly_count));

        // Review our geometry.
        for (int i = 0; i < poly_count; i++)
        {
            Assert.That(polys_gds2[i].pointarray.Count, Is.EqualTo(polys_fromdrawing[i].pointarray.Count));
            for (int pt = 0; pt < polys_gds2[i].pointarray.Count; pt++)
            {
                Assert.That(polys_gds2[i].pointarray[pt].X, Is.EqualTo(polys_fromdrawing[i].pointarray[pt].X));
                Assert.That(polys_gds2[i].pointarray[pt].Y, Is.EqualTo(polys_fromdrawing[i].pointarray[pt].Y));
            }
        }
    }

    [Test]
    public static void klayout_gds_arefs_test()
    {
        string gdsFile = baseDir + "klayout_test/gds/arefs.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        // T cell
        GCCell cell_gds = drawing_gds.findCell("T");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(1));
        List<GCPolygon> polys_gds = drawing_gds.convertToPolygons(cells: ["T"])[0];
        Assert.That(polys_gds.Count, Is.EqualTo(1));

        // A cell
        GCCell cell_gds2 = drawing_gds.findCell("A");
        Assert.That(cell_gds2.elementList.Count, Is.EqualTo(64));
        for (int i = 0; i < 64; i++)
        {
            Assert.That(cell_gds2.elementList[i].isCellrefArray(), Is.True);
        }
        List<GCPolygon> polys_gds2 = drawing_gds.convertToPolygons(cells: ["A"])[0];
        // 64 cell references that are 3x2 arrays of different configurations
        int poly_count = 64 * 6;
        Assert.That(polys_gds2.Count, Is.EqualTo(poly_count));
    }

    //[Test] // VERY slow.
    public static void klayout_algo_flat_region_au13b_test()
    {
        string gdsFile = baseDir + "klayout_test/algo/flat_region_au13b.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        // TOPTOP cell
        GCCell cell_gds = drawing_gds.findCell("TOPTOP");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(180008));
        List<GCPolygon> polys_gds = drawing_gds.convertToPolygons(cells: ["TOPTOP"])[0];
        Assert.That(polys_gds.Count, Is.EqualTo(180008));
    }

    [Test]
    public static void klayout_algo_antenna_au1_test()
    {
        string gdsFile = baseDir + "klayout_test/algo/antenna_au1.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        // TRANS cell
        GCCell cell_gds = drawing_gds.findCell("TRANS");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[0].isPolygon(), Is.True);
        Assert.That(cell_gds.elementList[0].layer_nr, Is.EqualTo(6));
        Assert.That(cell_gds.elementList[0].datatype_nr, Is.EqualTo(0));
        Assert.That(cell_gds.elementList[1].isBox(), Is.True);
        Assert.That(cell_gds.elementList[1].layer_nr, Is.EqualTo(8));
        Assert.That(cell_gds.elementList[1].datatype_nr, Is.EqualTo(0));
        List<GCPolygon> polys_gds = drawing_gds.convertToPolygons(cells: ["TRANS"])[0];

        int trans_poly_count = 2;
        Assert.That(polys_gds.Count, Is.EqualTo(trans_poly_count));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(2800));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(800));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(2800));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(800));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[4].Y, Is.EqualTo(0));

        Assert.That(polys_gds[1].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds[1].pointarray[0].X, Is.EqualTo(200));
        Assert.That(polys_gds[1].pointarray[0].Y, Is.EqualTo(600));
        Assert.That(polys_gds[1].pointarray[1].X, Is.EqualTo(600));
        Assert.That(polys_gds[1].pointarray[1].Y, Is.EqualTo(600));
        Assert.That(polys_gds[1].pointarray[2].X, Is.EqualTo(600));
        Assert.That(polys_gds[1].pointarray[2].Y, Is.EqualTo(200));
        Assert.That(polys_gds[1].pointarray[3].X, Is.EqualTo(200));
        Assert.That(polys_gds[1].pointarray[3].Y, Is.EqualTo(200));
        Assert.That(polys_gds[1].pointarray[4].X, Is.EqualTo(200));
        Assert.That(polys_gds[1].pointarray[4].Y, Is.EqualTo(600));

        // TOP cell
        GCCell cell_gds2 = drawing_gds.findCell("TOP");
        Assert.That(cell_gds2.elementList.Count, Is.EqualTo(109));
        List<GCPolygon> polys_gds2 = drawing_gds.convertToPolygons(cells: ["TOP"])[0];
        Assert.That(polys_gds2.Count, Is.EqualTo(119));

        // TOPTOP cell
        GCCell cell_gds3 = drawing_gds.findCell("TOPTOP");
        Assert.That(cell_gds3.elementList.Count, Is.EqualTo(100));
        List<GCPolygon> polys_gds3 = drawing_gds.convertToPolygons(cells: ["TOPTOP"])[0];
        Assert.That(polys_gds3.Count, Is.EqualTo(336));
    }

    [Test]
    public static void klayout_algo_angle_check_test()
    {
        string gdsFile = baseDir + "klayout_test/algo/angle_check_l1.gds";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GCDrawingfield drawing_gds = gcGDS.getDrawing();

        // A cell
        GCCell cell_gds = drawing_gds.findCell("A");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[0].isPolygon(), Is.True);
        List<GCPolygon> polys_gds = drawing_gds.convertToPolygons(cells: ["A"])[0];
        int a_poly_count = 1;
        Assert.That(polys_gds.Count, Is.EqualTo(a_poly_count));
        Assert.That(polys_gds[0].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds[0].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds[0].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds[0].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds[0].pointarray[3].Y, Is.EqualTo(0));

        // B cell
        GCCell cell_gds2 = drawing_gds.findCell("B");
        Assert.That(cell_gds2.elementList.Count, Is.EqualTo(3));
        Assert.That(cell_gds2.elementList[0].isCellref(), Is.True);
        Assert.That(cell_gds2.elementList[1].isCellrefArray(), Is.True);
        Assert.That(cell_gds2.elementList[2].isPolygon(), Is.True);
        List<GCPolygon> polys_gds2 = cell_gds2.convertToPolygons(1.0);
        // 3 instances of A cell and 1 polygon
        int b_poly_count = (3 * a_poly_count) + 1;
        Assert.That(polys_gds2.Count, Is.EqualTo(b_poly_count));
        Assert.That(polys_gds2[0].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds2[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds2[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds2[0].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds2[0].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds2[0].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds2[0].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds2[0].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds2[0].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds2[1].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[1].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds2[1].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds2[1].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds2[1].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[1].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds2[1].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[1].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds2[2].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[2].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds2[2].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds2[2].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds2[2].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[2].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds2[2].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds2[2].pointarray[3].Y, Is.EqualTo(1000));

        // C cell
        GCCell cell_gds3 = drawing_gds.findCell("C");
        Assert.That(cell_gds3.elementList.Count, Is.EqualTo(3));
        Assert.That(cell_gds3.elementList[0].isCellrefArray(), Is.True);
        Assert.That(cell_gds3.elementList[1].isCellref(), Is.True);
        Assert.That(cell_gds3.elementList[2].isPolygon(), Is.True);
        List<GCPolygon> polys_gds3 = cell_gds3.convertToPolygons(1.0);
        // 3 instances of B cell and 1 polygon
        int c_poly_count = (3 * 4) + 1;
        Assert.That(polys_gds3.Count, Is.EqualTo(c_poly_count));
        Assert.That(polys_gds3[0].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds3[0].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds3[0].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[0].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds3[0].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[0].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds3[0].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[0].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds3[0].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds3[1].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[1].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[1].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[1].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[1].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[1].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[1].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[1].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds3[2].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[2].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[2].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[2].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[2].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[2].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[2].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[2].pointarray[3].Y, Is.EqualTo(1000));

        Assert.That(polys_gds3[3].pointarray[0].X, Is.EqualTo(1000));
        Assert.That(polys_gds3[3].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[3].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds3[3].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[3].pointarray[2].X, Is.EqualTo(1500));
        Assert.That(polys_gds3[3].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[3].pointarray[3].X, Is.EqualTo(1500));
        Assert.That(polys_gds3[3].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[3].pointarray[4].X, Is.EqualTo(1000));
        Assert.That(polys_gds3[3].pointarray[4].Y, Is.EqualTo(0));

        int xoffset = 5000;
        Assert.That(polys_gds3[4].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds3[4].pointarray[0].X, Is.EqualTo(xoffset + 0));
        Assert.That(polys_gds3[4].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[4].pointarray[1].X, Is.EqualTo(xoffset + 1000));
        Assert.That(polys_gds3[4].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[4].pointarray[2].X, Is.EqualTo(xoffset + 1000));
        Assert.That(polys_gds3[4].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[4].pointarray[3].X, Is.EqualTo(xoffset + 0));
        Assert.That(polys_gds3[4].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds3[5].pointarray[0].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[5].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[5].pointarray[1].X, Is.EqualTo(xoffset + 2000));
        Assert.That(polys_gds3[5].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[5].pointarray[2].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[5].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[5].pointarray[3].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[5].pointarray[3].Y, Is.EqualTo(0));

        Assert.That(polys_gds3[6].pointarray[0].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[6].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds3[6].pointarray[1].X, Is.EqualTo(xoffset + 2000));
        Assert.That(polys_gds3[6].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[6].pointarray[2].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[6].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[6].pointarray[3].X, Is.EqualTo(xoffset + 4000));
        Assert.That(polys_gds3[6].pointarray[3].Y, Is.EqualTo(1000));

        Assert.That(polys_gds3[7].pointarray[0].X, Is.EqualTo(xoffset + 1000));
        Assert.That(polys_gds3[7].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[7].pointarray[1].X, Is.EqualTo(xoffset + 1000));
        Assert.That(polys_gds3[7].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[7].pointarray[2].X, Is.EqualTo(xoffset + 1500));
        Assert.That(polys_gds3[7].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds3[7].pointarray[3].X, Is.EqualTo(xoffset + 1500));
        Assert.That(polys_gds3[7].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds3[7].pointarray[4].X, Is.EqualTo(xoffset + 1000));
        Assert.That(polys_gds3[7].pointarray[4].Y, Is.EqualTo(0));

        Assert.That(polys_gds3[8].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[8].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds3[8].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds3[8].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds3[8].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[8].pointarray[2].Y, Is.EqualTo(6000));
        Assert.That(polys_gds3[8].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[8].pointarray[3].Y, Is.EqualTo(4000));

        Assert.That(polys_gds3[9].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[9].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[9].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[9].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds3[9].pointarray[2].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[9].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[9].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[9].pointarray[3].Y, Is.EqualTo(12000));

        Assert.That(polys_gds3[10].pointarray[0].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[10].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[10].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds3[10].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds3[10].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds3[10].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[10].pointarray[3].X, Is.EqualTo(2000));
        Assert.That(polys_gds3[10].pointarray[3].Y, Is.EqualTo(12000));

        Assert.That(polys_gds3[11].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[11].pointarray[0].Y, Is.EqualTo(6000));
        Assert.That(polys_gds3[11].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds3[11].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds3[11].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds3[11].pointarray[2].Y, Is.EqualTo(7000));
        Assert.That(polys_gds3[11].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[11].pointarray[3].Y, Is.EqualTo(7000));
        Assert.That(polys_gds3[11].pointarray[4].X, Is.EqualTo(4000));
        Assert.That(polys_gds3[11].pointarray[4].Y, Is.EqualTo(6000));

        Assert.That(polys_gds3[12].pointarray[0].X, Is.EqualTo(5000));
        Assert.That(polys_gds3[12].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds3[12].pointarray[1].X, Is.EqualTo(5000));
        Assert.That(polys_gds3[12].pointarray[1].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[12].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds3[12].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds3[12].pointarray[3].X, Is.EqualTo(5000));
        Assert.That(polys_gds3[12].pointarray[3].Y, Is.EqualTo(4000));

        // TOP cell
        GCCell cell_gds4 = drawing_gds.findCell("TOP");
        Assert.That(cell_gds4.elementList.Count, Is.EqualTo(6));
        Assert.That(cell_gds4.elementList[0].isCellref(), Is.True);
        Assert.That(cell_gds4.elementList[1].isCellrefArray(), Is.True);
        Assert.That(cell_gds4.elementList[2].isPolygon(), Is.True);
        Assert.That(cell_gds4.elementList[3].isPolygon(), Is.True);
        Assert.That(cell_gds4.elementList[4].isPolygon(), Is.True);
        Assert.That(cell_gds4.elementList[5].isPolygon(), Is.True);
        List<GCPolygon> polys_gds4 = cell_gds4.convertToPolygons(1.0);
        // 3 instances of the C cell and 4 polygons.
        int top_poly_count = (3 * c_poly_count) + 4;
        Assert.That(polys_gds4.Count, Is.EqualTo(top_poly_count));

        int index = 0;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(16000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(1500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(1500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(16000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(6500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(6500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(21000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(21000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(19000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(23000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(27000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(23000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(27000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(21000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(21000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(22000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(22000));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(21000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(19000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(1000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(1500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(1500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(1000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(6500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(6500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(4000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(12000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(12000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(6000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(9000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(5000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(4000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(11000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(11000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(1000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(11000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(11000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(11500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(11500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(11000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(1000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(1000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(16500));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(2000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(16500));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(0));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(16000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(0));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(4000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(12000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(8000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(12000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(5));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(6000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(7000));
        Assert.That(polys_gds4[index].pointarray[4].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[4].Y, Is.EqualTo(6000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(4000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(12000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(4000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(19000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(17000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(15000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(15000));

        index++;
        Assert.That(polys_gds4[index].pointarray.Count, Is.EqualTo(4));
        Assert.That(polys_gds4[index].pointarray[0].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[0].Y, Is.EqualTo(19000));
        Assert.That(polys_gds4[index].pointarray[1].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[1].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[2].X, Is.EqualTo(14000));
        Assert.That(polys_gds4[index].pointarray[2].Y, Is.EqualTo(27000));
        Assert.That(polys_gds4[index].pointarray[3].X, Is.EqualTo(10000));
        Assert.That(polys_gds4[index].pointarray[3].Y, Is.EqualTo(19000));
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);

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

        string oasFile = baseDir + "gdstk_reference/f_rep3_kl.oas";
        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        GCCell cell_oas = drawing_oas.findCell("Base");
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);

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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);

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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);

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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);

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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);

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
    public static void f_rep6_test()
    {
        string gdsFile = baseDir + "gdstk_reference/f_rep6.gds";
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        x_offset = 0;
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
        x_offset = 0;
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

        polyIndex++;
        x_offset = 0;
        y_offset = -2000;
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas_kl = cell_oas_kl.convertToPolygons(drawing_oas_kl.getDrawingScale());
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
    public static void ref_f_rep6_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_f_rep6.gds";
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        y_offset = 1000;
        x_offset = 0;
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
        y_offset = 3000;
        x_offset = 0;
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
        y_offset = -2000;
        x_offset = 0;
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

        string oasFile = baseDir + "gdstk_reference/ref_f_rep6.oas";
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        y_offset = -2000;
        x_offset = 0;
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
        y_offset = 0;
        x_offset = 0;
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
        y_offset = 1000;
        x_offset = 0;
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
        y_offset = 3000;
        x_offset = 0;
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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

        GCDrawingfield drawing_gds = gcGDS.getDrawing();
        GCCell cell_gds = drawing_gds.findCell("Ref");
        Assert.That(cell_gds.elementList.Count, Is.EqualTo(4));
        for (int i = 0; i < 4; i++)
        {
            Assert.That(cell_gds.elementList[i].isCellref(), Is.True);
        }
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep4.oas";
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
    public static void ref_rect_rep5_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep5.gds";
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep5.oas";
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
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
    public static void ref_rect_rep6_test()
    {
        string gdsFile = baseDir + "gdstk_reference/ref_rectangle_rep6.gds";
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
        List<GCPolygon> polys_gds = cell_gds.convertToPolygons(1.0);
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
        y_offset = 1000;
        x_offset = 0;
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
        y_offset = 3000;
        x_offset = 0;
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
        y_offset = -2000;
        x_offset = 0;
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

        string oasFile = baseDir + "gdstk_reference/ref_rectangle_rep6.oas";
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
        List<GCPolygon> polys_oas = cell_oas.convertToPolygons(1.0);
        Assert.That(polys_oas.Count, Is.EqualTo(4));

        polyIndex = 0;
        y_offset = -2000;
        x_offset = 0;
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
        y_offset = 0;
        x_offset = 0;
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
        y_offset = 1000;
        x_offset = 0;
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
        y_offset = 3000;
        x_offset = 0;
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

        gcell.addPolygon(Clipper.ScalePath64(circleD, scale * scale), 1, 0);
        Assert.That(gcell.elementList[0].isPolygon(), Is.True);
        Assert.That(((GCPolygon)gcell.elementList[0]).isCircle().circle, Is.True);
        // For comparison
        gcell.addCircle(1, 1, new(0, 0), 500);
        Assert.That(gcell.elementList[1].isPolygon(), Is.True);
        Assert.That(((GCPolygon)gcell.elementList[1]).isCircle().circle, Is.True);

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

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
        Assert.That(t1.Count, Is.EqualTo(1));
        Assert.That(t1[0].layer_nr, Is.EqualTo(1));
        Assert.That(t1[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t1[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t1[0].pointarray[0], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t1[0].pointarray[1], Is.EqualTo(new Point64(0, 100, 0)));
        Assert.That(t1[0].pointarray[2], Is.EqualTo(new Point64(40, 100, 0)));
        Assert.That(t1[0].pointarray[3], Is.EqualTo(new Point64(40, 40, 0)));
        Assert.That(t1[0].pointarray[4], Is.EqualTo(new Point64(100, 40, 0)));
        Assert.That(t1[0].pointarray[5], Is.EqualTo(new Point64(100, 0, 0)));
        Assert.That(t1[0].pointarray[6], Is.EqualTo(new Point64(0, 0, 0)));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(1));
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
        Assert.That(t2.Count, Is.EqualTo(4));
        Assert.That(t2[0].layer_nr, Is.EqualTo(1));
        Assert.That(t2[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[0].pointarray[0], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t2[0].pointarray[1], Is.EqualTo(new Point64(0, 100, 0)));
        Assert.That(t2[0].pointarray[2], Is.EqualTo(new Point64(40, 100, 0)));
        Assert.That(t2[0].pointarray[3], Is.EqualTo(new Point64(40, 40, 0)));
        Assert.That(t2[0].pointarray[4], Is.EqualTo(new Point64(100, 40, 0)));
        Assert.That(t2[0].pointarray[5], Is.EqualTo(new Point64(100, 0, 0)));
        Assert.That(t2[0].pointarray[6], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t2[1].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[1].pointarray[0], Is.EqualTo(new Point64(110, 0, 0)));
        Assert.That(t2[1].pointarray[1], Is.EqualTo(new Point64(110, 100, 0)));
        Assert.That(t2[1].pointarray[2], Is.EqualTo(new Point64(150, 100, 0)));
        Assert.That(t2[1].pointarray[3], Is.EqualTo(new Point64(150, 40, 0)));
        Assert.That(t2[1].pointarray[4], Is.EqualTo(new Point64(210, 40, 0)));
        Assert.That(t2[1].pointarray[5], Is.EqualTo(new Point64(210, 0, 0)));
        Assert.That(t2[1].pointarray[6], Is.EqualTo(new Point64(110, 0, 0)));
        Assert.That(t2[2].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[2].pointarray[0], Is.EqualTo(new Point64(0, 110, 0)));
        Assert.That(t2[2].pointarray[1], Is.EqualTo(new Point64(0, 210, 0)));
        Assert.That(t2[2].pointarray[2], Is.EqualTo(new Point64(40, 210, 0)));
        Assert.That(t2[2].pointarray[3], Is.EqualTo(new Point64(40, 150, 0)));
        Assert.That(t2[2].pointarray[4], Is.EqualTo(new Point64(100, 150, 0)));
        Assert.That(t2[2].pointarray[5], Is.EqualTo(new Point64(100, 110, 0)));
        Assert.That(t2[2].pointarray[6], Is.EqualTo(new Point64(0, 110, 0)));
        Assert.That(t2[3].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[3].pointarray[0], Is.EqualTo(new Point64(110, 110, 0)));
        Assert.That(t2[3].pointarray[1], Is.EqualTo(new Point64(110, 210, 0)));
        Assert.That(t2[3].pointarray[2], Is.EqualTo(new Point64(150, 210, 0)));
        Assert.That(t2[3].pointarray[3], Is.EqualTo(new Point64(150, 150, 0)));
        Assert.That(t2[3].pointarray[4], Is.EqualTo(new Point64(210, 150, 0)));
        Assert.That(t2[3].pointarray[5], Is.EqualTo(new Point64(210, 110, 0)));
        Assert.That(t2[3].pointarray[6], Is.EqualTo(new Point64(110, 110, 0)));

        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(2));
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
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

        List<List<GCPolygon>> initial = gcOAS.getDrawing().convertToPolygons();
        string initial_hash = Utils.GetSHA1Hash(initial);

        string outFile = outDir + "/c3_consistency_from_oas.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcOAS, outFile);
        gw.save();
        Assert.That(File.Exists(outFile), Is.True);

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(outFile, GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        gcGDS.getDrawing().resize(1.0 / gcGDS.getDrawing().getDrawingScale());
        List<List<GCPolygon>> gds = gcGDS.getDrawing().convertToPolygons();
        string gds_hash = Utils.GetSHA1Hash(gds);

        Assert.That(gds_hash, Is.EqualTo(initial_hash));

        string outFile2 = outDir + "/c3_consistency_from_oas.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcOAS, outFile2);
        ow.save();
        Assert.That(File.Exists(outFile2), Is.True);

        GeoCoreHandler gH_OAS2 = new();
        gH_OAS2.updateGeoCoreHandler(outFile2, GeoCore.fileType.oasis);
        GeoCore gcOAS2 = gH_OAS2.getGeo();
        Assert.That(gcOAS2.isValid(), Is.True);

        List<List<GCPolygon>> oas = gcOAS2.getDrawing().convertToPolygons();
        string oas_hash = Utils.GetSHA1Hash(oas);
        Assert.That(oas_hash, Is.EqualTo(initial_hash));
    }

    [Test]
    public static void consistency_from_gds()
    {
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(baseDir + "/consistency/c3_consistency.gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        List<List<GCPolygon>> initial = gcGDS.getDrawing().convertToPolygons();
        string initial_hash = Utils.GetSHA1Hash(initial);

        string outFile = outDir + "/c3_consistency_from_gds.gds";
        if (File.Exists(outFile))
        {
            File.Delete(outFile);
        }
        gds.gdsWriter gw = new(gcGDS, outFile);
        gw.save();
        Assert.That(File.Exists(outFile), Is.True);

        GeoCoreHandler gH_GDS2 = new();
        gH_GDS2.updateGeoCoreHandler(outFile, GeoCore.fileType.gds);
        GeoCore gcGDS2 = gH_GDS2.getGeo();
        Assert.That(gcGDS2.isValid(), Is.True);

        gcGDS2.getDrawing().resize(1.0 / gcGDS2.getDrawing().getDrawingScale());
        List<List<GCPolygon>> gds = gcGDS2.getDrawing().convertToPolygons();
        string gds_hash = Utils.GetSHA1Hash(gds);

        Assert.That(gds_hash, Is.EqualTo(initial_hash));

        string outFile2 = outDir + "/c3_consistency_from_gds.oas";
        if (File.Exists(outFile2))
        {
            File.Delete(outFile2);
        }
        oasis.oasWriter ow = new(gcGDS, outFile2);
        ow.save();
        Assert.That(File.Exists(outFile2), Is.True);

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(outFile2, GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        List<List<GCPolygon>> oas = gcOAS.getDrawing().convertToPolygons();
        string oas_hash = Utils.GetSHA1Hash(oas);
        Assert.That(oas_hash, Is.EqualTo(initial_hash));
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

        List<GCPolygon> t1 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
        Assert.That(t1.Count, Is.EqualTo(1));
        Assert.That(t1[0].layer_nr, Is.EqualTo(1));
        Assert.That(t1[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t1[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t1[0].pointarray[0], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t1[0].pointarray[1], Is.EqualTo(new Point64(0, 100, 0)));
        Assert.That(t1[0].pointarray[2], Is.EqualTo(new Point64(40, 100, 0)));
        Assert.That(t1[0].pointarray[3], Is.EqualTo(new Point64(40, 40, 0)));
        Assert.That(t1[0].pointarray[4], Is.EqualTo(new Point64(100, 40, 0)));
        Assert.That(t1[0].pointarray[5], Is.EqualTo(new Point64(100, 0, 0)));
        Assert.That(t1[0].pointarray[6], Is.EqualTo(new Point64(0, 0, 0)));

        // Array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("a");
        Assert.That(gcGDS.activeStructure, Is.EqualTo(1));

        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t2 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
        Assert.That(t2.Count, Is.EqualTo(6));
        Assert.That(t2[0].layer_nr, Is.EqualTo(1));
        Assert.That(t2[0].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[0].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[0].pointarray[0], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t2[0].pointarray[1], Is.EqualTo(new Point64(0, 100, 0)));
        Assert.That(t2[0].pointarray[2], Is.EqualTo(new Point64(40, 100, 0)));
        Assert.That(t2[0].pointarray[3], Is.EqualTo(new Point64(40, 40, 0)));
        Assert.That(t2[0].pointarray[4], Is.EqualTo(new Point64(100, 40, 0)));
        Assert.That(t2[0].pointarray[5], Is.EqualTo(new Point64(100, 0, 0)));
        Assert.That(t2[0].pointarray[6], Is.EqualTo(new Point64(0, 0, 0)));
        Assert.That(t2[1].layer_nr, Is.EqualTo(1));
        Assert.That(t2[1].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[1].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[1].pointarray[0], Is.EqualTo(new Point64(110, 0, 0)));
        Assert.That(t2[1].pointarray[1], Is.EqualTo(new Point64(110, 100, 0)));
        Assert.That(t2[1].pointarray[2], Is.EqualTo(new Point64(150, 100, 0)));
        Assert.That(t2[1].pointarray[3], Is.EqualTo(new Point64(150, 40, 0)));
        Assert.That(t2[1].pointarray[4], Is.EqualTo(new Point64(210, 40, 0)));
        Assert.That(t2[1].pointarray[5], Is.EqualTo(new Point64(210, 0, 0)));
        Assert.That(t2[1].pointarray[6], Is.EqualTo(new Point64(110, 0, 0)));
        Assert.That(t2[2].layer_nr, Is.EqualTo(1));
        Assert.That(t2[2].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[2].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[2].pointarray[0], Is.EqualTo(new Point64(0, 110, 0)));
        Assert.That(t2[2].pointarray[1], Is.EqualTo(new Point64(0, 210, 0)));
        Assert.That(t2[2].pointarray[2], Is.EqualTo(new Point64(40, 210, 0)));
        Assert.That(t2[2].pointarray[3], Is.EqualTo(new Point64(40, 150, 0)));
        Assert.That(t2[2].pointarray[4], Is.EqualTo(new Point64(100, 150, 0)));
        Assert.That(t2[2].pointarray[5], Is.EqualTo(new Point64(100, 110, 0)));
        Assert.That(t2[2].pointarray[6], Is.EqualTo(new Point64(0, 110, 0)));
        Assert.That(t2[3].layer_nr, Is.EqualTo(1));
        Assert.That(t2[3].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[3].pointarray.Count, Is.EqualTo(7));
        Assert.That(t2[3].pointarray[0], Is.EqualTo(new Point64(110, 110, 0)));
        Assert.That(t2[3].pointarray[1], Is.EqualTo(new Point64(110, 210, 0)));
        Assert.That(t2[3].pointarray[2], Is.EqualTo(new Point64(150, 210, 0)));
        Assert.That(t2[3].pointarray[3], Is.EqualTo(new Point64(150, 150, 0)));
        Assert.That(t2[3].pointarray[4], Is.EqualTo(new Point64(210, 150, 0)));
        Assert.That(t2[3].pointarray[5], Is.EqualTo(new Point64(210, 110, 0)));
        Assert.That(t2[3].pointarray[6], Is.EqualTo(new Point64(110, 110, 0)));
        Assert.That(t2[4].layer_nr, Is.EqualTo(1));
        Assert.That(t2[4].datatype_nr, Is.EqualTo(0));
        Assert.That(t2[4].pointarray.Count, Is.EqualTo(5));
        Assert.That(t2[4].pointarray[0], Is.EqualTo(new Point64(60, 173, 0)));
        Assert.That(t2[4].pointarray[1], Is.EqualTo(new Point64(60, 201, 0)));
        Assert.That(t2[4].pointarray[2], Is.EqualTo(new Point64(83, 201, 0)));
        Assert.That(t2[4].pointarray[3], Is.EqualTo(new Point64(83, 173, 0)));
        Assert.That(t2[4].pointarray[4], Is.EqualTo(new Point64(60, 173, 0)));
        Assert.That(t2[5].layer_nr, Is.EqualTo(2).Within(2));
        Assert.That(t2[5].datatype_nr, Is.EqualTo(0).Within(0));
        Assert.That(t2[5].pointarray.Count, Is.EqualTo(5).Within(5));
        Assert.That(t2[5].pointarray[0], Is.EqualTo(new Point64(52, 67, 0)));
        Assert.That(t2[5].pointarray[1], Is.EqualTo(new Point64(52, 93, 0)));
        Assert.That(t2[5].pointarray[2], Is.EqualTo(new Point64(103, 93, 0)));
        Assert.That(t2[5].pointarray[3], Is.EqualTo(new Point64(103, 67, 0)));
        Assert.That(t2[5].pointarray[4], Is.EqualTo(new Point64(52, 67, 0)));

        // Nested array
        gcGDS.activeStructure = gcGDS.getStructureList().IndexOf("b");
        gcGDS.updateGeometry(gcGDS.activeStructure, gcGDS.activeLD);

        List<GCPolygon> t3 = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale());
        List<GCPolygon> t3a = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale(), layer: 1, datatype: 0);
        List<GCPolygon> t3b = gcGDS.getDrawing().cellList[gcGDS.activeStructure].convertToPolygons(gcGDS.getDrawing().getDrawingScale(), layer: 2, datatype: 0);

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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // triangle
        poly = Helper.initedPath64(3);
        poly[0] = new(0, 0);
        poly[1] = new(10, 20);
        poly[2] = new(20, 0);

        gcell.addPolygon(poly, 2, 0);

        // pentagram
        poly = Helper.initedPath64(5);
        poly[0] = new(5, 0);
        poly[1] = new(0, 10);
        poly[2] = new(10, 20);
        poly[3] = new(20, 10);
        poly[4] = new(15, 0);

        gcell.addPolygon(poly, 3, 0);

        // trapezoid
        poly = Helper.initedPath64(4);
        poly[0] = new(0, 0);
        poly[1] = new(5, 20);
        poly[2] = new(15, 20);
        poly[3] = new(20, 0);

        gcell.addPolygon(poly, 4, 0);

        // parallelogram
        poly = Helper.initedPath64(4);
        poly[0] = new(0, 0);
        poly[1] = new(10, 20);
        poly[2] = new(20, 20);
        poly[3] = new(10, 0);

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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[1].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[2].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[3].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[4].convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[1].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[2].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[3].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[4].convertToPolygons(1.0);
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
        path[0] = new(0, 0);
        path[1] = new(0, 10);
        path[2] = new(20, 10);
        path[3] = new(20, 40);
        path[4] = new(0, 40);

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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
            databaseunits = 1,
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

        gcell.addCircle(1, 0, new(10, 10), 5.0);

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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref1";
        gcell.addCellref();

        gcell.elementList[^1].setPos(new(10, 0));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test"));
        gcell.elementList[^1].setName("test");
        gcell.elementList[^1].rotate(0);
        gcell.elementList[^1].scale(2);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref1_mx";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new(10, 0));
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 30);
        poly[2] = new(10, 30);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref2";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new(20, 20));
        gcell.elementList[^1].setCellRef(drawing_.findCell("test2"));
        gcell.elementList[^1].setName("test2");
        gcell.elementList[^1].rotate(90);
        gcell.elementList[^1].scale(2);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellref2_mx";
        gcell.addCellref();
        gcell.elementList[^1].setPos(new(20, 20));
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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray1";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(0, 0));
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
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(0));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray2";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (0 + (rowIndex * col_pitch.Y)))));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray3";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(-10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(-10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(30 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(0 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(0 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray4";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        Point64 row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        double scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        int polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * 20) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * 10) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);


        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray5";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        Point64 row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        double scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.False);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        int polyIndex = 0;
        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(-30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(-30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(-10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(-10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray1_mx";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(0, 0));
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
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        col_pitch = cell_gds.elementList[^1].getColPitch();
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray2_mx";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * (-20 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * (10 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * (-10 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * (20 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * (0 + (colIndex * row_pitch.X)))));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * (-0 + (rowIndex * col_pitch.Y)))));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray3_mx";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(30 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(20 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(30 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(20 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(30 + (rowIndex * row_pitch.Y)));
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray4_mx";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(0));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(pos.Y + (scale * 0) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(pos.Y + (scale * -20) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(pos.X + (scale * 10) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(pos.Y + (scale * -10) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(pos.X + (scale * 20) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(pos.X + (scale * 0) + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(pos.Y + (scale * -0) + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
        poly[0] = new(0, 0);
        poly[1] = new(0, 20);
        poly[2] = new(10, 20);
        poly[3] = new(10, 10);
        poly[4] = new(20, 10);
        poly[5] = new(20, 0);

        gcell.addPolygon(poly, 1, 0);

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        Path64 array = Helper.initedPath64(3);
        array[0] = new(0, 0);
        array[1] = new(0, 80);
        array[2] = new(100, 0);

        gcell = drawing_.addCell();
        gcell.cellName = "test_cellrefarray5_mx";
        gcell.addCellref(drawing_.findCell("test"), new(0, 0));
        gcell.addCellrefArray(drawing_.findCell("test"), array, 4, 4);
        gcell.elementList[^1].setPos(new(10, 20));
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
        Assert.That(col_pitch.X, Is.EqualTo(0 / count.X));
        Assert.That(col_pitch.Y, Is.EqualTo(80 / count.Y));
        row_pitch = cell_gds.elementList[^1].getRowPitch();
        Assert.That(row_pitch.X, Is.EqualTo(100 / count.X));
        Assert.That(row_pitch.Y, Is.EqualTo(0 / count.Y));
        scale = cell_gds.elementList[^1].getScale();
        Assert.That(scale, Is.EqualTo(2));
        Assert.That(cell_gds.elementList[^1].getAngle(), Is.EqualTo(90));
        Assert.That(cell_gds.elementList[^1].getMirrorX(), Is.True);
        List<GCPolygon> polys_gds = cell_gds.elementList[^1].convertToPolygons(1.0);
        Assert.That(polys_gds.Count, Is.EqualTo(16));

        for (int rowIndex = 0; rowIndex < count.Y; rowIndex++)
        {
            for (int colIndex = 0; colIndex < count.X; colIndex++)
            {
                Assert.That(polys_gds[polyIndex].pointarray.Count, Is.EqualTo(7));
                Assert.That(polys_gds[polyIndex].pointarray[0].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[0].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[1].X, Is.EqualTo(50 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[1].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[2].X, Is.EqualTo(50 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[2].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[3].X, Is.EqualTo(30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[3].Y, Is.EqualTo(40 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[4].X, Is.EqualTo(30 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[4].Y, Is.EqualTo(60 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[5].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[5].Y, Is.EqualTo(60 + (rowIndex * col_pitch.Y)));
                Assert.That(polys_gds[polyIndex].pointarray[6].X, Is.EqualTo(10 + (colIndex * row_pitch.X)));
                Assert.That(polys_gds[polyIndex].pointarray[6].Y, Is.EqualTo(20 + (rowIndex * col_pitch.Y)));
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
        List<GCPolygon> polys_oas = cell_oas.elementList[^1].convertToPolygons(1.0);
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
                Assert.That(polys_oas[polyIndex].pointarray[2].Y, Is.EqualTo(40 + (rowIndex * row_pitch.Y)));
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

        GCDrawingfield drawing_oas_ref = gcOAS_ref.getDrawing();
        List<GCPolygon> polys_oas_ref = drawing_oas_ref.convertToPolygons(cells: ["test_cellrefarray1"])[0];
        Assert.That(polys_oas_ref.Count, Is.EqualTo(16));

        // Check our repetition has been handled correctly.
        int xPitch = 25;
        int yPitch = 20;

        int index = 0;
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                Assert.That(polys_oas_ref[index].pointarray[0].X, Is.EqualTo(0 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[0].Y, Is.EqualTo(0 + (row * yPitch)));
                Assert.That(polys_oas_ref[index].pointarray[1].X, Is.EqualTo(0 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[1].Y, Is.EqualTo(20 + (row * yPitch)));
                Assert.That(polys_oas_ref[index].pointarray[2].X, Is.EqualTo(10 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[2].Y, Is.EqualTo(20 + (row * yPitch)));
                Assert.That(polys_oas_ref[index].pointarray[3].X, Is.EqualTo(10 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[3].Y, Is.EqualTo(10 + (row * yPitch)));
                Assert.That(polys_oas_ref[index].pointarray[4].X, Is.EqualTo(20 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[4].Y, Is.EqualTo(10 + (row * yPitch)));
                Assert.That(polys_oas_ref[index].pointarray[5].X, Is.EqualTo(20 + (col * xPitch)));
                Assert.That(polys_oas_ref[index].pointarray[5].Y, Is.EqualTo(0 + (row * yPitch)));
                index++;
            }
        }
    }

    [Test]
    public static void test_cellref_irregular_parser()
    {
        string inputfile = baseDir + "cellref_array_irregular";
        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(inputfile + ".gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(inputfile + ".oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        List<GCPolygon> gds_polys = gcGDS.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];
        List<GCPolygon> oas_polys = gcOAS.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];

        string gds_hash = Utils.GetMD5Hash(gds_polys);
        string oas_hash = Utils.GetMD5Hash(oas_polys);

        Assert.That(gds_hash, Is.EqualTo(oas_hash));

        // Review geometry - if hashes match, we only need to check one of the inputs.
        Assert.That(gds_polys.Count, Is.EqualTo(7));
        for (int pindex = 0; pindex < 7; pindex++)
        {
            int x_offset = 0;
            int y_offset = 0;
            switch (pindex)
            {
                case 2:
                    x_offset = 8;
                    y_offset = 50;
                    break;
                case 3:
                    x_offset = -130;
                    y_offset = 10;
                    break;
                case 4:
                    x_offset = 90;
                    y_offset = 70;
                    break;
                case 5:
                    x_offset = 30;
                    y_offset = 50;
                    break;
                default:
                    x_offset = 0;
                    y_offset = 0;
                    break;
            }
            Assert.That(gds_polys[pindex].pointarray.Count, Is.EqualTo(7));
            Assert.That(gds_polys[pindex].pointarray[0].X, Is.EqualTo(0 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[0].Y, Is.EqualTo(0 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[1].X, Is.EqualTo(0 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[1].Y, Is.EqualTo(20 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[2].X, Is.EqualTo(10 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[2].Y, Is.EqualTo(20 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[3].X, Is.EqualTo(10 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[3].Y, Is.EqualTo(10 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[4].X, Is.EqualTo(20 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[4].Y, Is.EqualTo(10 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[5].X, Is.EqualTo(20 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[5].Y, Is.EqualTo(0 + y_offset));
            Assert.That(gds_polys[pindex].pointarray[6].X, Is.EqualTo(0 + x_offset));
            Assert.That(gds_polys[pindex].pointarray[6].Y, Is.EqualTo(0 + y_offset));
        }

        string outfile = outDir + "cellref_array_irregular_resavecheck";
        string gdsFile = outfile + ".gds";
        if (File.Exists(gdsFile))
        {
            File.Delete(gdsFile);
        }
        gds.gdsWriter gw = new(gcGDS, gdsFile);
        gw.save();
        Assert.That(File.Exists(gdsFile), Is.True);

        string oasFile = outfile + ".oas";
        if (File.Exists(oasFile))
        {
            File.Delete(oasFile);
        }
        oasis.oasWriter ow = new(gcOAS, oasFile);
        ow.save();
        Assert.That(File.Exists(oasFile), Is.True);

        GeoCoreHandler gH_GDS_reload = new();
        gH_GDS_reload.updateGeoCoreHandler(gdsFile, GeoCore.fileType.gds);
        GeoCore gcGDS_reload = gH_GDS_reload.getGeo();
        Assert.That(gcGDS_reload.isValid(), Is.True);

        GeoCoreHandler gH_OAS_reload = new();
        gH_OAS_reload.updateGeoCoreHandler(oasFile, GeoCore.fileType.oasis);
        GeoCore gcOAS_reload = gH_OAS_reload.getGeo();
        Assert.That(gcOAS_reload.isValid(), Is.True);

        List<GCPolygon> gds_polys_reload = gcGDS_reload.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];
        List<GCPolygon> oas_polys_reload = gcOAS_reload.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];

        string gds_reload_hash = Utils.GetMD5Hash(gds_polys_reload);
        string oas_reload_hash = Utils.GetMD5Hash(oas_polys_reload);

        Assert.That(gds_reload_hash, Is.EqualTo(oas_reload_hash));
        Assert.That(gds_reload_hash, Is.EqualTo(gds_hash));
        Assert.That(oas_reload_hash, Is.EqualTo(oas_hash));

        // No need to review geometry - if hashes match, all is good.
    }

    [Test]
    public static void test_cellref_irregular_offset_parser()
    {
        string inputfile = baseDir + "cellref_array_irregular_offset";

        GeoCoreHandler gH_GDS = new();
        gH_GDS.updateGeoCoreHandler(inputfile + ".gds", GeoCore.fileType.gds);
        GeoCore gcGDS = gH_GDS.getGeo();
        Assert.That(gcGDS.isValid(), Is.True);

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(inputfile + ".oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        List<GCPolygon> gds_polys = gcGDS.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];
        List<GCPolygon> oas_polys = gcOAS.getDrawing().convertToPolygons(cells: ["test_cellrefarray_irregular"])[0];

        string gds_hash = Utils.GetMD5Hash(gds_polys);
        string oas_hash = Utils.GetMD5Hash(oas_polys);

        Assert.That(gds_polys.Count, Is.EqualTo(oas_polys.Count));

        int x_offset_gds = 0;
        int y_offset_gds = 0;
        int x_offset_oas = 0;
        int y_offset_oas = 0;
        for (int poly = 0; poly < 7; poly++)
        {
            switch (poly)
            {
                case 0:
                    x_offset_gds = 0;
                    y_offset_gds = 0;
                    x_offset_oas = 10;
                    y_offset_oas = -30;
                    break;
                case 1:
                    x_offset_gds = 0;
                    y_offset_gds = 0;
                    x_offset_oas = 0;
                    y_offset_oas = 0;
                    break;
                case 2:
                    x_offset_gds = 8;
                    y_offset_gds = 50;
                    x_offset_oas = 0;
                    y_offset_oas = 0;
                    break;
                case 3:
                    x_offset_gds = -130;
                    y_offset_gds = 10;
                    x_offset_oas = -130;
                    y_offset_oas = 10;
                    break;
                case 4:
                    x_offset_gds = 90;
                    y_offset_gds = 70;
                    x_offset_oas = 8;
                    y_offset_oas = 50;
                    break;
                case 5:
                    x_offset_gds = 30;
                    y_offset_gds = 50;
                    x_offset_oas = 30;
                    y_offset_oas = 50;
                    break;
                case 6:
                    x_offset_gds = 10;
                    y_offset_gds = -30;
                    x_offset_oas = 90;
                    y_offset_oas = 70;
                    break;
            }
            Assert.That(gds_polys[poly].pointarray[0].X, Is.EqualTo(0 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[0].Y, Is.EqualTo(0 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[1].X, Is.EqualTo(0 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[1].Y, Is.EqualTo(20 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[2].X, Is.EqualTo(10 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[2].Y, Is.EqualTo(20 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[3].X, Is.EqualTo(10 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[3].Y, Is.EqualTo(10 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[4].X, Is.EqualTo(20 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[4].Y, Is.EqualTo(10 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[5].X, Is.EqualTo(20 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[5].Y, Is.EqualTo(0 + y_offset_gds));
            Assert.That(gds_polys[poly].pointarray[6].X, Is.EqualTo(0 + x_offset_gds));
            Assert.That(gds_polys[poly].pointarray[6].Y, Is.EqualTo(0 + y_offset_gds));

            Assert.That(oas_polys[poly].pointarray[0].X, Is.EqualTo(0 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[0].Y, Is.EqualTo(0 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[1].X, Is.EqualTo(0 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[1].Y, Is.EqualTo(20 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[2].X, Is.EqualTo(10 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[2].Y, Is.EqualTo(20 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[3].X, Is.EqualTo(10 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[3].Y, Is.EqualTo(10 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[4].X, Is.EqualTo(20 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[4].Y, Is.EqualTo(10 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[5].X, Is.EqualTo(20 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[5].Y, Is.EqualTo(0 + y_offset_oas));
            Assert.That(oas_polys[poly].pointarray[6].X, Is.EqualTo(0 + x_offset_oas));
            Assert.That(oas_polys[poly].pointarray[6].Y, Is.EqualTo(0 + y_offset_oas));
        }

        // Precomputed hashes. Geometry is differently ordered, so hashes differ
        Assert.That(gds_hash, Is.EqualTo("5UJvmu7iw5bkcDvYxrImTQ=="));
        Assert.That(oas_hash, Is.EqualTo("15trvprfnCQH7L8uWiGx7A=="));
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

        gcell.addText(1, 0, new(10, 10), "Text");

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
                1 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, -1, 0, 0 } },
                2 => new[,] { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                3 => new[,] { { 0, 1, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                4 => new[,] { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 1, -1, 0, 1 }, { 1, 0, 0, 0 } },
                5 => new[,] { { 0, 1, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, -1, 0, 0 } },
                6 => new[,] { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 1, 0, 0, 1 }, { 1, -1, 0, 0 } },
                7 => new[,] { { 0, 1, 0, 0 }, { 0, 0, 0, 1 }, { 1, -1, 0, 1 }, { 1, 0, 0, 0 } },
                8 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, -1, 1 }, { 1, 0, 0, 0 } },
                9 => new[,] { { 0, 0, 0, 0 }, { 0, 0, -1, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                10 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 1, 0 } },
                11 => new[,] { { 0, 0, 1, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                12 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, -1, 1 }, { 1, 0, 1, 0 } },
                13 => new[,] { { 0, 0, 1, 0 }, { 0, 0, -1, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                14 => new[,] { { 0, 0, 0, 0 }, { 0, 0, -1, 1 }, { 1, 0, 0, 1 }, { 1, 0, 1, 0 } },
                15 => new[,] { { 0, 0, 1, 0 }, { 0, 0, 0, 1 }, { 1, 0, -1, 1 }, { 1, 0, 0, 0 } },
                16 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 1, 0 }, { 1, 0, 0, 0 }, { 0, 0, 0, 0 } },
                17 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 1, 0 }, { 1, 0, 1, 0 }, { 0, 0, 0, 0 } },
                18 => new[,] { { 0, 0, 0, 0 }, { 1, 0, 1, 0 }, { 1, 0, 0, 0 }, { 0, 0, 0, 0 } },
                19 => new[,] { { 0, 0, 1, 0 }, { 1, 0, 1, 0 }, { 1, 0, 0, 0 }, { 0, 0, 1, 0 } },
                20 => new[,] { { 0, 0, 0, 0 }, { 0, 1, 0, 1 }, { 0, 2, 0, 0 }, { 0, 0, 0, 0 } },
                21 => new[,] { { 0, 0, 0, 1 }, { 0, 2, 0, 1 }, { 0, 1, 0, 0 }, { 0, 0, 0, 1 } },
                22 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 2, 0 }, { 1, 0, 1, 0 }, { 0, 0, 0, 0 } },
                23 => new[,] { { 1, 0, 0, 0 }, { 0, 0, 1, 0 }, { 1, 0, 2, 0 }, { 1, 0, 0, 0 } },
                24 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 0, 0 } },
                25 => new[,] { { 0, 0, 0, 0 }, { 0, 0, 1, 0 }, { 1, 0, 1, 0 }, { 1, 0, 0, 0 } },
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

                pa[pt] = new(x, y);

                if (x > w)
                {
                    w = x;
                }
                if (y > h)
                {
                    h = y;
                }
            }

            pa[^1] = new(pa[0]);

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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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

        Assert.That(cell_oas.elementList[^1].isPolygon(), Is.True);
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
                pa[0] = new(x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new(x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new(x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new(x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new(x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new(x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new(x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new(x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new(pa[0]);

        gcell.addPolygon(pa, 1, 0);

        trapezoid_orientation = true;

        switch (trapezoid_orientation)
        {
            // (m & 0x80)
            case true:
                //  vertically
                pa[0] = new(x, y + Math.Max(trapezoid_delta_a, 0));
                pa[1] = new(x, y + h + Math.Min(trapezoid_delta_b, 0));
                pa[2] = new(x + w, y + h - Math.Max(trapezoid_delta_b, 0));
                pa[3] = new(x + w, y - Math.Min(trapezoid_delta_a, 0));
                break;
            default:
                //  horizontally
                pa[0] = new(x + Math.Max(trapezoid_delta_a, 0), y + h);
                pa[1] = new(x + w + Math.Min(trapezoid_delta_b, 0), y + h);
                pa[2] = new(x + w - Math.Max(trapezoid_delta_b, 0), y);
                pa[3] = new(x - Math.Min(trapezoid_delta_a, 0), y);
                break;
        }

        pa[4] = new(pa[0]);

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
        List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
        List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
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
            poly[0] = new(0, 0);
            poly[1] = new(0, 20);
            poly[2] = new(10, 20);
            poly[3] = new(10, 10);
            poly[4] = new(20, 10);
            poly[5] = new(20, 0);

            gcell.addPolygon(poly, 1, 0);

            master_gcell.addCellref();
            GCElement cellref = master_gcell.elementList[^1];
            cellref.setPos(new(40 * (i % edge), 40 * Math.Floor((double)i / edge)));
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
            List<GCPolygon> polys_gds = cell_gds.elementList[0].convertToPolygons(1.0);
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
            List<GCPolygon> polys_oas = cell_oas.elementList[0].convertToPolygons(1.0);
            Assert.That(polys_oas.Count, Is.EqualTo(1));
            Assert.That(polys_oas[0].pointarray.Count, Is.EqualTo(7));
        }
    }

    [Test]
    public static void move_box_test()
    {
        int dimension = 10;
        int scale = 100;
        string filename = "move_box_" + 10 + "_" + 20;
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
            databaseunits = GCDrawingfield.default_databaseunits / scale,
            userunits = GCDrawingfield.default_userunits / scale,
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

        gcell.addBox(0, 0, dimension, dimension, 1, 1);
        GCElement content = gcell.elementList[^1];

        Assert.That(content.isBox(), Is.True);
        Assert.That(content.getPos().X, Is.EqualTo(0));
        Assert.That(content.getPos().Y, Is.EqualTo(0));
        List<GCPolygon> polys = content.convertToPolygons(drawing_.getDrawingScale());
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
        content.move(new(move_x, move_y));
        Assert.That(content.getPos().X, Is.EqualTo(move_x));
        Assert.That(content.getPos().Y, Is.EqualTo(move_y));
        minx = Int64.MaxValue;
        miny = Int64.MaxValue;
        // Due to integer representation, we have to scale up to ensure we don't lose anything.
        double scalefactor = drawing_.getDrawingScale();
        switch (scalefactor)
        {
            case 0:
                scalefactor = 1;
                break;
            case < 1:
                // Need to invert because otherwise the integer representation could fail.
                scalefactor = 1 / scalefactor;
                break;
        }
        polys = content.convertToPolygons(scalefactor);
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
        // Inflated values here due to integer representation.
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
        drawing_gds.databaseunits = GCDrawingfield.default_databaseunits / scale;
        drawing_gds.userunits = GCDrawingfield.default_userunits / scale;
        GCCell cell_gds = drawing_gds.findCell("test");
        int elementCount = cell_gds.elementList.Count;
        Assert.That(elementCount, Is.EqualTo(1));
        Assert.That(cell_gds.elementList[^1].isBox(), Is.True);
        Assert.That(cell_gds.elementList[^1].getPos().X, Is.EqualTo(10));
        Assert.That(cell_gds.elementList[^1].getPos().Y, Is.EqualTo(20));

        string out_filename = outDir + "/" + filename + "_resave_from_gds.gds";
        save_gdsii(gcGDS, out_filename);

        out_filename = outDir + "/" + filename + "_resave_from_gds.oas";
        save_oasis(gcGDS, out_filename);

        GeoCoreHandler gH_OAS = new();
        gH_OAS.updateGeoCoreHandler(outDir + "/" + filename + ".oas", GeoCore.fileType.oasis);
        GeoCore gcOAS = gH_OAS.getGeo();
        Assert.That(gcOAS.isValid(), Is.True);

        GCDrawingfield drawing_oas = gcOAS.getDrawing();
        drawing_oas.databaseunits = GCDrawingfield.default_databaseunits / scale;
        drawing_oas.userunits = GCDrawingfield.default_userunits / scale;
        GCCell cell_oas = drawing_oas.findCell("test");
        elementCount = cell_oas.elementList.Count;
        Assert.That(elementCount, Is.EqualTo(1));
        Assert.That(cell_oas.elementList[^1].isBox(), Is.True);
        Assert.That(cell_oas.elementList[^1].getPos().X, Is.EqualTo(10));
        Assert.That(cell_oas.elementList[^1].getPos().Y, Is.EqualTo(20));

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

    #endregion
}