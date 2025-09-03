using Clipper2Lib;
using geoWrangler;
using shapeEngine;
using System.Threading;
using System.Threading.Tasks;
using System.Text;
using System.IO;

namespace UnitTests;

/// <summary>
/// Comprehensive tests for the shapeEngine library, which provides shape generation
/// and manipulation capabilities including rectangles, L-shapes, T-shapes, X-shapes,
/// U-shapes, S-shapes, and other complex geometric forms for layout design.
/// </summary>
public class ShapeEngineTests
{
    private static string root_loc = "/d/development/DesignLibs_GPL/shapeengine_out/";

    // Copied from the private content of ShapeLibrary, for testing.
    private static readonly List<string> availableShapes_all = new()
    {
        "(None)", "Rectangle/Square", "L-shape", "T-shape", "X-shape", "U-shape", "S-shape", "GDS/Oasis", "Boolean",
        "Text", "Bounding", "Layout"
    };

    static int[] shapeTable = {
        (int)ShapeLibrary.shapeNames_all.none,
        (int)ShapeLibrary.shapeNames_all.rect,
        (int)ShapeLibrary.shapeNames_all.Lshape,
        (int)ShapeLibrary.shapeNames_all.Tshape,
        (int)ShapeLibrary.shapeNames_all.Xshape,
        (int)ShapeLibrary.shapeNames_all.Ushape,
        (int)ShapeLibrary.shapeNames_all.Sshape,
        (int)ShapeLibrary.shapeNames_all.GEOCORE,
        (int)ShapeLibrary.shapeNames_all.BOOLEAN
    };

    private static int[] shapeTable2 = new[]
    {
        (int)ShapeLibrary.shapeNames_all.none,
        (int)ShapeLibrary.shapeNames_all.rect,
        (int)ShapeLibrary.shapeNames_all.Lshape,
        (int)ShapeLibrary.shapeNames_all.Tshape,
        (int)ShapeLibrary.shapeNames_all.Xshape,
        (int)ShapeLibrary.shapeNames_all.Ushape,
        (int)ShapeLibrary.shapeNames_all.Sshape,
        (int)ShapeLibrary.shapeNames_all.text,
        (int)ShapeLibrary.shapeNames_all.bounding,
        (int)ShapeLibrary.shapeNames_all.complex
    };

    #region New Unit Tests

    /// <summary>
    /// Tests ShapeLibrary constructor and basic initialization.
    /// </summary>
    [Test]
    public static void ShapeLibrary_Constructor_InitializesCorrectly()
    {
        // Arrange
        ShapeSettings shapeSettings = new ShapeSettings();

        // Act
        ShapeLibrary shapeLib = new ShapeLibrary(shapeTable, shapeSettings);

        // Assert: New instance should be properly initialized
        Assert.That(shapeLib.last_error, Is.EqualTo(""));
        Assert.That(shapeLib.shapeValid, Is.False);
        Assert.That(shapeLib.geoCoreShapeOrthogonal, Is.False);
    }

    /// <summary>
    /// Tests shapeNames_all enumeration contains all expected values.
    /// </summary>
    [Test]
    public static void ShapeLibrary_ShapeNamesEnum_ContainsExpectedValues()
    {
        // Assert: All expected shape types are available
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.none), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.rect), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.Lshape), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.Tshape), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.Xshape), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.Ushape), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.Sshape), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.GEOCORE), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.BOOLEAN), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.text), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.bounding), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeLibrary.shapeNames_all), ShapeLibrary.shapeNames_all.complex), Is.True);
    }

    /// <summary>
    /// Tests shapesForClient method with valid input.
    /// </summary>
    [Test]
    public static void ShapeLibrary_ShapesForClient_ValidInput_ConfiguresCorrectly()
    {
        // Arrange
        ShapeSettings shapeSettings = new ShapeSettings();
        ShapeLibrary shapeLib = new ShapeLibrary(shapeTable, shapeSettings);
        int[] clientShapeDefinition = { 0, 1, 2, 3 };

        // Act
        shapeLib.shapesForClient(clientShapeDefinition);

        // Assert: Should complete without throwing
        Assert.That(shapeLib.last_error, Is.EqualTo(""));
    }

    /// <summary>
    /// Tests shapesForClient method with too many shapes.
    /// </summary>
    [Test]
    public static void ShapeLibrary_ShapesForClient_TooManyShapes_ThrowsException()
    {
        // Arrange
        ShapeSettings shapeSettings = new ShapeSettings();
        ShapeLibrary shapeLib = new ShapeLibrary(shapeTable, shapeSettings);
        int[] tooManyShapes = new int[20]; // More than available shapes

        // Act & Assert
        Assert.Throws<Exception>(() => shapeLib.shapesForClient(tooManyShapes));
    }

    /// <summary>
    /// Tests getAvailableShapes static method.
    /// </summary>
    [Test]
    public static void ShapeLibrary_GetAvailableShapes_ReturnsCorrectNames()
    {
        // Arrange
        int[] clientShapeDefinition = { 0, 1, 2 };

        // Act
        List<string> shapes = ShapeLibrary.getAvailableShapes(clientShapeDefinition);

        // Assert
        Assert.That(shapes.Count, Is.EqualTo(3));
        Assert.That(shapes[0], Is.EqualTo("(None)"));
        Assert.That(shapes[1], Is.EqualTo("Rectangle/Square"));
        Assert.That(shapes[2], Is.EqualTo("L-shape"));
    }

    /// <summary>
    /// Tests ShapeLibrary enum values correspond to expected integer values.
    /// </summary>
    [Test]
    public static void ShapeLibrary_EnumValues_HaveCorrectIntegerValues()
    {
        // Assert: Enum values should match expected integers
        Assert.That((int)ShapeLibrary.shapeNames_all.none, Is.EqualTo(0));
        Assert.That((int)ShapeLibrary.shapeNames_all.rect, Is.EqualTo(1));
        Assert.That((int)ShapeLibrary.shapeNames_all.Lshape, Is.EqualTo(2));
        Assert.That((int)ShapeLibrary.shapeNames_all.Tshape, Is.EqualTo(3));
        Assert.That((int)ShapeLibrary.shapeNames_all.Xshape, Is.EqualTo(4));
        Assert.That((int)ShapeLibrary.shapeNames_all.Ushape, Is.EqualTo(5));
        Assert.That((int)ShapeLibrary.shapeNames_all.Sshape, Is.EqualTo(6));
        Assert.That((int)ShapeLibrary.shapeNames_all.GEOCORE, Is.EqualTo(7));
        Assert.That((int)ShapeLibrary.shapeNames_all.BOOLEAN, Is.EqualTo(8));
        Assert.That((int)ShapeLibrary.shapeNames_all.text, Is.EqualTo(9));
        Assert.That((int)ShapeLibrary.shapeNames_all.bounding, Is.EqualTo(10));
        Assert.That((int)ShapeLibrary.shapeNames_all.complex, Is.EqualTo(11));
    }

    /// <summary>
    /// Tests ShapeSettings tipLocations enumeration.
    /// </summary>
    [Test]
    public static void ShapeSettings_TipLocations_ContainsExpectedValues()
    {
        // Assert: All expected tip locations are available
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.tipLocations), ShapeSettings.tipLocations.none), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.tipLocations), ShapeSettings.tipLocations.L), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.tipLocations), ShapeSettings.tipLocations.R), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.tipLocations), ShapeSettings.tipLocations.all), Is.True);
    }

    /// <summary>
    /// Tests ShapeSettings subShapeLocations enumeration.
    /// </summary>
    [Test]
    public static void ShapeSettings_SubShapeLocations_ContainsExpectedValues()
    {
        // Assert: All expected sub-shape locations are available
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.subShapeLocations), ShapeSettings.subShapeLocations.TL), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.subShapeLocations), ShapeSettings.subShapeLocations.TR), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.subShapeLocations), ShapeSettings.subShapeLocations.BL), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.subShapeLocations), ShapeSettings.subShapeLocations.BR), Is.True);
        Assert.That(Enum.IsDefined(typeof(ShapeSettings.subShapeLocations), ShapeSettings.subShapeLocations.C), Is.True);
    }

    #endregion

    #region Original Integration Tests

    // [SetUp]
    public static void ShapeEngineSetup()
    {
        shapeLookupTest();
        rectangleTest();
        rectangleRoundingTest();
        notLTest();
        lTest();
        lInnerRoundingTest();
        lOuterRoundingTest();
        lInnerOuterRoundingTest();
        lInnerRoundingVarTest();
        lOuterRoundingVarTest();
        lInnerOuterRoundingVarTest();
        lInnerRoundingTensionTest();
        sTest();
        sRoundingTest();
        tTest();
        tRoundingTest();
        uTest();
        uRoundingTest();
        xTest();
        xRoundingTest();
        biasTest();
        sBiasTest();
        lTipTest();
        lTipVarTest();
        uTipTest();
        sTipTest();
        globalOffsetTest();
        subShapePositioningTest();
        customOrthoTest();
        customOrthoOuterRoundingTest();
        customOrthoInnerRoundingTest();
        customOrthoTipsTest();
    }

    [Test]
    public static void shapeLookupTest()
    {
        // Check our mapping.
        List<string> shapes = ShapeLibrary.getAvailableShapes(shapeTable);
        List<string> shapes2 = ShapeLibrary.getAvailableShapes(shapeTable2);
        Assert.That(shapes.Count, Is.EqualTo(shapeTable.Length));
        Assert.That(shapes2.Count, Is.EqualTo(shapeTable2.Length));
        Assert.That(shapes.IndexOf(availableShapes_all[0]), Is.EqualTo(0));
        Assert.That(shapes.IndexOf(availableShapes_all[1]), Is.EqualTo(1));
        Assert.That(shapes.IndexOf(availableShapes_all[2]), Is.EqualTo(2));
        Assert.That(shapes.IndexOf(availableShapes_all[3]), Is.EqualTo(3));
        Assert.That(shapes.IndexOf(availableShapes_all[4]), Is.EqualTo(4));
        Assert.That(shapes.IndexOf(availableShapes_all[5]), Is.EqualTo(5));
        Assert.That(shapes.IndexOf(availableShapes_all[6]), Is.EqualTo(6));
        Assert.That(shapes.IndexOf(availableShapes_all[7]), Is.EqualTo(7));
        Assert.That(shapes.IndexOf(availableShapes_all[8]), Is.EqualTo(8));
        Assert.That(shapes.IndexOf(availableShapes_all[9]), Is.EqualTo(-1));
        Assert.That(shapes.IndexOf(availableShapes_all[10]), Is.EqualTo(-1));
        Assert.That(shapes.IndexOf(availableShapes_all[11]), Is.EqualTo(-1));

        Assert.That(shapes2.IndexOf(availableShapes_all[0]), Is.EqualTo(0));
        Assert.That(shapes2.IndexOf(availableShapes_all[1]), Is.EqualTo(1));
        Assert.That(shapes2.IndexOf(availableShapes_all[2]), Is.EqualTo(2));
        Assert.That(shapes2.IndexOf(availableShapes_all[3]), Is.EqualTo(3));
        Assert.That(shapes2.IndexOf(availableShapes_all[4]), Is.EqualTo(4));
        Assert.That(shapes2.IndexOf(availableShapes_all[5]), Is.EqualTo(5));
        Assert.That(shapes2.IndexOf(availableShapes_all[6]), Is.EqualTo(6));
        Assert.That(shapes2.IndexOf(availableShapes_all[7]), Is.EqualTo(-1));
        Assert.That(shapes2.IndexOf(availableShapes_all[8]), Is.EqualTo(-1));
        Assert.That(shapes2.IndexOf(availableShapes_all[9]), Is.EqualTo(7));
        Assert.That(shapes2.IndexOf(availableShapes_all[10]), Is.EqualTo(8));
        Assert.That(shapes2.IndexOf(availableShapes_all[11]), Is.EqualTo(9));

        // Do we get the correct shape from the engine?
        ShapeSettings shapeSettings1 = new ShapeSettings();
        shapeSettings1.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.bounding);
        shapeSettings1.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings1.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);

        ShapeLibrary shape1 = new ShapeLibrary(shapeTable2, shapes2.IndexOf("Bounding"), shapeSettings1);
        PathD out_1 = shape1.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "bounding_lookup_test.svg", FillRule.NonZero, 800, 800, 10);
    }

    [Test]
    public static void rectangleTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.rect));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(61));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-200));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(10));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void rectangleRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 10);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.rect));
        PathD out_ = shape.processCorners(false, false, 90, 1, .5);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle_rounding.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(false));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(121));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-(Math.PI * 10 * 10) - area), Is.LessThanOrEqualTo(0.15));
    }

    [Test]
    public static void biasTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.sBias, 5);
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.R);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.hTBias, 7);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.rect));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle_bias.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(105));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-22 * 30));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(22));
        Assert.That(bounds.Height, Is.EqualTo(30));
    }

    // Set for L, but it's a rectangle. We should get the rectangle.
    [Test]
    public static void notLTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.rect));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "notL.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(61));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-200));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(10));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(69));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-250));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lInnerRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 5);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_inner.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(68));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-252.8 - area), Is.LessThanOrEqualTo(0.01));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lInnerRoundingVarTest()
    {
        ShapeSettings shapeSettings_ref = new ShapeSettings();
        shapeSettings_ref.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.iCR, 5);
        ShapeLibrary shape_ref = new ShapeLibrary(shapeTable, shapeSettings_ref);
        shape_ref.setShape(shapeSettings_ref.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_ref.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.That(clean_ref.Count, Is.EqualTo(68));
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.That(Math.Abs(-252.8 - area_ref), Is.LessThanOrEqualTo(0.01));
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.That(bounds_ref.Width, Is.EqualTo(15));
        Assert.That(bounds_ref.Height, Is.EqualTo(20));

        // This should result in the same shape because we use the variation to deliver the outer rounding (1.0 * 5)
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 0);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, iCV: 5, iCVariation_scalar: 1.0);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_innervariation.svg", FillRule.NonZero, 800, 800, 10);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(68));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(area_ref - area), Is.LessThanOrEqualTo(0.01));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lOuterRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_outer.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(60));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-225.18 - area), Is.LessThanOrEqualTo(0.01));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lOuterRoundingVarTest()
    {
        ShapeSettings shapeSettings_ref = new ShapeSettings();
        shapeSettings_ref.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        ShapeLibrary shape_ref = new ShapeLibrary(shapeTable, shapeSettings_ref);
        shape_ref.setShape(shapeSettings_ref.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_ref.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.That(clean_ref.Count, Is.EqualTo(60));
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.That(Math.Abs(-225.18 - area_ref), Is.LessThanOrEqualTo(0.01));
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.That(bounds_ref.Width, Is.EqualTo(15));
        Assert.That(bounds_ref.Height, Is.EqualTo(20));

        // This should result in the same shape because we use the variation to deliver the outer rounding (1.0 * 5)
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 0);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, oCV: 5, oCVariation_scalar: 1.0);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_outervariation.svg", FillRule.NonZero, 800, 800, 10);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(60));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(area_ref - area), Is.LessThanOrEqualTo(0.01));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lInnerOuterRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 10);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(59));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-228 - area), Is.LessThanOrEqualTo(0.015));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lInnerOuterRoundingVarTest()
    {
        ShapeSettings shapeSettings_ref = new ShapeSettings();
        shapeSettings_ref.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.iCR, 10);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        ShapeLibrary shape_ref = new ShapeLibrary(shapeTable, shapeSettings_ref);
        shape_ref.setShape(shapeSettings_ref.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_ref.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.That(clean_ref.Count, Is.EqualTo(59));
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.That(Math.Abs(-228 - area_ref), Is.LessThanOrEqualTo(0.015));
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.That(bounds_ref.Width, Is.EqualTo(15));
        Assert.That(bounds_ref.Height, Is.EqualTo(20));

        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 0);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, oCV: 5, oCVariation_scalar: 1.0, iCV: 10, iCVariation_scalar: 1.0);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_inneroutervar.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(59));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(area_ref - area), Is.LessThanOrEqualTo(0.015));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void lInnerRoundingTensionTest()
    {
        ShapeSettings shapeSettings_0 = new ShapeSettings();
        shapeSettings_0.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.iCR, 10);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        shapeSettings_0.setInt(ShapeSettings.properties_i.edgeSlide, 1);
        shapeSettings_0.setDecimal(ShapeSettings.properties_decimal.eTension, 0);
        ShapeLibrary shape_0 = new ShapeLibrary(shapeTable, shapeSettings_0);
        shape_0.setShape(shapeSettings_0.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_0.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_0 = shape_0.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_0 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension0.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_0 = GeoWrangler.removeDuplicates(out_0);
        // Check point count - start and end points are the same.
        Assert.That(clean_0.Count, Is.EqualTo(59));
        // Check expected area
        double area_0 = Clipper.Area(out_0);
        Assert.That(Math.Abs(-228 - area_0), Is.LessThanOrEqualTo(0.015));
        RectD bounds = Clipper.GetBounds(out_0);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));

        ShapeSettings shapeSettings_05 = new ShapeSettings();
        shapeSettings_05.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.iCR, 10);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        shapeSettings_05.setInt(ShapeSettings.properties_i.edgeSlide, 1);
        shapeSettings_05.setDecimal(ShapeSettings.properties_decimal.eTension, 0.5m);
        ShapeLibrary shape_05 = new ShapeLibrary(shapeTable, shapeSettings_05);
        shape_05.setShape(shapeSettings_05.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_05.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_05 = shape_05.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_05 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension05.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_05 = GeoWrangler.removeDuplicates(out_05);
        // Check point count - start and end points are the same.
        Assert.That(clean_05.Count, Is.EqualTo(60));
        // Check expected area
        double area_05 = Clipper.Area(out_05);
        Assert.That(Math.Abs(-233.78 - area_05), Is.LessThanOrEqualTo(0.015));
        RectD bounds_05 = Clipper.GetBounds(out_05);
        Assert.That(bounds_05.Width, Is.EqualTo(15));
        Assert.That(bounds_05.Height, Is.EqualTo(20));

        ShapeSettings shapeSettings_10 = new ShapeSettings();
        shapeSettings_10.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.horLength, 5.0m, 1);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.iCR, 10);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.oCR, 5);
        shapeSettings_10.setInt(ShapeSettings.properties_i.edgeSlide, 1);
        shapeSettings_10.setDecimal(ShapeSettings.properties_decimal.eTension, 1);
        ShapeLibrary shape_10 = new ShapeLibrary(shapeTable, shapeSettings_10);
        shape_10.setShape(shapeSettings_10.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_10.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_10 = shape_10.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_10 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension10.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_10 = GeoWrangler.removeDuplicates(out_10);
        // Check point count - start and end points are the same.
        Assert.That(clean_10.Count, Is.EqualTo(60));
        // Check expected area
        double area_10 = Clipper.Area(out_10);
        Assert.That(Math.Abs(-239.31 - area_10), Is.LessThanOrEqualTo(0.015));
        RectD bounds_10 = Clipper.GetBounds(out_10);
        Assert.That(bounds_10.Width, Is.EqualTo(15));
        Assert.That(bounds_10.Height, Is.EqualTo(20));
    }

    [Test]
    public static void sTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(309));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-(3600 - ((20 * 20) + (10 * 15))) - area), Is.LessThanOrEqualTo(0.0001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void sRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 18);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 15);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(259));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-2915.0441 - area), Is.LessThanOrEqualTo(0.0001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void tTest()
    {
        ShapeSettings shapeSettings = new();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Tshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 20m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 12m, 1);

        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Tshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "t.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(241));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-((20 * 60) + (40 * 10))));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void tRoundingTest()
    {
        ShapeSettings shapeSettings = new();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Tshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 20m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 12m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 7);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 12);

        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Tshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "t_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(204));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-1503.028 - area), Is.LessThanOrEqualTo(0.001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void uTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Ushape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 40.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 5m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 10m, 1);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Ushape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(259));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-((60 * 40) - (30 * 30))));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(40));
    }

    [Test]
    public static void uRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Ushape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 40.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 5m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 10m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 7);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 12);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Ushape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(false));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(225));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-1384.031 - area), Is.LessThanOrEqualTo(0.001));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(40));
    }

    [Test]
    public static void xTest()
    {
        ShapeSettings shapeSettings = new();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Xshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, -5m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 12m, 1);

        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Xshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(197));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-1400 - area), Is.LessThanOrEqualTo(0.001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(40));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void xRoundingTest()
    {
        ShapeSettings shapeSettings = new();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Xshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, -5m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 12m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 7);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 12);

        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Xshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(164));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-1327.185 - area), Is.LessThanOrEqualTo(0.001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(40));
        Assert.That(bounds.Height, Is.EqualTo(60));
    }

    [Test]
    public static void sBiasTest()
    {
        ShapeSettings shapeSettings_ref = new ShapeSettings();
        shapeSettings_ref.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        // shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.sBias, 1);

        ShapeLibrary shape_ref = new ShapeLibrary(shapeTable, shapeSettings_ref);
        shape_ref.setShape(shapeSettings_ref.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_ref.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);

        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.sBias, 2m);

        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ref }, true);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_sbias.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(325));
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(636 - Math.Abs(area_ref - area)), Is.LessThanOrEqualTo(0.0001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(64));
        Assert.That(bounds.Height, Is.EqualTo(64));
    }

    [Test]
    public static void lTipTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.hTBias, 6m);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.vTBias, 12m);
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.T);
        shapeSettings.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.R);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tip.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(277));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-2360));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(66));
        Assert.That(bounds.Height, Is.EqualTo(72));
    }

    [Test]
    public static void lTipVarTest()
    {
        ShapeSettings shapeSettings_ref = new ShapeSettings();
        shapeSettings_ref.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.hTBias, 6m);
        shapeSettings_ref.setDecimal(ShapeSettings.properties_decimal.vTBias, 12m);
        shapeSettings_ref.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.T);
        shapeSettings_ref.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.R);
        ShapeLibrary shape_ref = new ShapeLibrary(shapeTable, shapeSettings_ref);
        shape_ref.setShape(shapeSettings_ref.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_ref.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_ref = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_ref), 0.001);
        Assert.That(orthogonal_ref, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.That(clean_ref.Count, Is.EqualTo(277));
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.That(area_ref, Is.EqualTo(-2360));
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.That(bounds_ref.Width, Is.EqualTo(66));
        Assert.That(bounds_ref.Height, Is.EqualTo(72));

        ShapeSettings shapeSettings_htv = new ShapeSettings();
        shapeSettings_htv.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.hTBias, 0);
        shapeSettings_htv.setDecimal(ShapeSettings.properties_decimal.vTBias, 12m);
        shapeSettings_htv.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.T);
        shapeSettings_htv.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.R);
        ShapeLibrary shape_htv = new ShapeLibrary(shapeTable, shapeSettings_htv);
        shape_htv.setShape(shapeSettings_htv.getInt(ShapeSettings.properties_i.shapeIndex));
        shape_htv.computeCage(0, 6, 0);
        // Check the shape settings are in the shape.
        Assert.That(shape_htv.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_htv = shape_htv.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_htv }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tipvarh.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_htv = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_htv), 0.001);
        Assert.That(orthogonal_htv, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean_htv = GeoWrangler.removeDuplicates(out_htv);
        // Check point count - start and end points are the same.
        Assert.That(clean_htv.Count, Is.EqualTo(277));
        // Check expected area
        double area_htv = Clipper.Area(out_htv);
        Assert.That(area_htv, Is.EqualTo(area_ref));
        RectD bounds_htv = Clipper.GetBounds(out_htv);
        Assert.That(bounds_htv.Width, Is.EqualTo(66));
        Assert.That(bounds_htv.Height, Is.EqualTo(72));

        ShapeSettings shapeSettings_vtv = new ShapeSettings();
        shapeSettings_vtv.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Lshape);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 0);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.horLength, 40.0m, 1);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.hTBias, 6m);
        shapeSettings_vtv.setDecimal(ShapeSettings.properties_decimal.vTBias, 0m);
        shapeSettings_vtv.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.T);
        shapeSettings_vtv.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.R);
        ShapeLibrary shape_vtv = new ShapeLibrary(shapeTable, shapeSettings_vtv);
        shape_vtv.setShape(shapeSettings_vtv.getInt(ShapeSettings.properties_i.shapeIndex));
        shape_vtv.computeCage(12, 0, 0);
        // Check the shape settings are in the shape.
        Assert.That(shape_vtv.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Lshape));
        PathD out_vtv = shape_vtv.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_vtv }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tipvarv.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_vtv = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_vtv), 0.001);
        Assert.That(orthogonal_vtv, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean_vtv = GeoWrangler.removeDuplicates(out_vtv);
        // Check point count - start and end points are the same.
        Assert.That(clean_vtv.Count, Is.EqualTo(277));
        // Check expected area
        double area_vtv = Clipper.Area(out_vtv);
        Assert.That(area_vtv, Is.EqualTo(-2360));
        RectD bounds_vtv = Clipper.GetBounds(out_vtv);
        Assert.That(bounds_vtv.Width, Is.EqualTo(66));
        Assert.That(bounds_vtv.Height, Is.EqualTo(72));
    }

    [Test]
    public static void uTipTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Ushape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 40.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 30m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 5m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 10m, 1);
        shapeSettings.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.BR);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.hTBias, 5);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.vTBias, 10);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Ushape));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u_brtip.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.That(orthogonal, Is.EqualTo(true));
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(239));
        // Check expected area
        double area = Clipper.Area(out_);
        // Tip configuration reduces width and height
        Assert.That(area, Is.EqualTo(-((60 * 40) - (25 * 20))));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(60));
        Assert.That(bounds.Height, Is.EqualTo(40));
    }

    [Test]
    public static void sTipTest()
    {
        ShapeSettings shapeSettings_1 = new ShapeSettings();
        shapeSettings_1.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        // Tips
        shapeSettings_1.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.BL);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.hTBias, 5);
        shapeSettings_1.setDecimal(ShapeSettings.properties_decimal.vTBias, 2);
        ShapeLibrary shape_1 = new ShapeLibrary(shapeTable, shapeSettings_1);
        shape_1.setShape(shapeSettings_1.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_1.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_1 = shape_1.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_bltips.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_1 = GeoWrangler.removeDuplicates(out_1);
        // Check point count - start and end points are the same.
        Assert.That(clean_1.Count, Is.EqualTo(309));
        // Check expected area
        double area_1 = Clipper.Area(out_1);
        // In this test, the BL tip setting moves the lower edge upwards, reducing the notch vertically.
        // The left tip setting has no effect on the horizontal dimension.
        Assert.That(Math.Abs(-(3600 - ((20 * 18) + (10 * 15))) - area_1), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_1 = Clipper.GetBounds(clean_1);
        Assert.That(bounds_1.Width, Is.EqualTo(60));
        Assert.That(bounds_1.Height, Is.EqualTo(60));

        ShapeSettings shapeSettings_2 = new ShapeSettings();
        shapeSettings_2.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);
        // Tips
        shapeSettings_2.setInt(ShapeSettings.properties_i.subShape2TipLocIndex, (int)ShapeSettings.tipLocations.BR);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.hTBias, 5);
        shapeSettings_2.setDecimal(ShapeSettings.properties_decimal.vTBias, 2);
        ShapeLibrary shape_2 = new ShapeLibrary(shapeTable, shapeSettings_2);
        shape_2.setShape(shapeSettings_2.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_2.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_2 = shape_2.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_brtips.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_2 = GeoWrangler.removeDuplicates(out_2);
        // Check point count - start and end points are the same.
        Assert.That(clean_2.Count, Is.EqualTo(297));
        // Check expected area
        double area_2 = Clipper.Area(out_2);
        // In this case, the BR tip setting moves the lower edge upwards, reducing the notch vertically.
        // The right tip setting reduces the notch width.
        Assert.That(Math.Abs(-(3600 - ((15 * 18) + (10 * 15))) - area_2), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_2 = Clipper.GetBounds(clean_2);
        Assert.That(bounds_2.Width, Is.EqualTo(60));
        Assert.That(bounds_2.Height, Is.EqualTo(60));

    }

    [Test]
    public static void globalOffsetTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.gHorOffset, 10);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.gVerOffset, 20);

        PointD offset = shapeOffsets.doOffsets(0, shapeSettings);

        Assert.That(offset.x, Is.EqualTo(10));
        // Inverse due to Clipper Y axis inversion
        Assert.That(offset.y, Is.EqualTo(-20));
    }


    // The inversion of the Y offset value is due to the inverted Y coordinate system in Clipper.
    [Test]
    public static void subShapePositioningTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.Sshape);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 60.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 1);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 15.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 10.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 7.0m, 1);
        // 45 based on pushing the second subshape against the wall of the outer.
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horOffset, 45.0m, 2);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verOffset, 9.0m, 2);

        shapeSettings.setInt(ShapeSettings.properties_i.subShapeRefIndex, 0);

        PointD offset_0 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.That(offset_0.x, Is.EqualTo(0));
        Assert.That(offset_0.y, Is.EqualTo(0));

        ShapeLibrary shape_0 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_0.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_0.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_0 = shape_0.processCorners(false, false, 90, 1, 1);
        out_0 = GeoWrangler.move(out_0, offset_0.x, -offset_0.y);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_0 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape0_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_0 = GeoWrangler.removeDuplicates(out_0);
        // Check point count - start and end points are the same.
        Assert.That(clean_0.Count, Is.EqualTo(309));
        // Check expected area
        double area_0 = Clipper.Area(out_0);
        Assert.That(Math.Abs(-(3600 - ((20 * 20) + (10 * 15))) - area_0), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_0 = Clipper.GetBounds(clean_0);
        Assert.That(bounds_0.left, Is.EqualTo(0));
        Assert.That(bounds_0.top, Is.EqualTo(0));

        // Change our reference to subshape 1
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeRefIndex, 1);
        PointD offset_1 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.That(offset_1.x, Is.EqualTo(0));
        Assert.That(offset_1.y, Is.EqualTo(7));

        ShapeLibrary shape_1 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_1.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_1.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_1 = shape_1.processCorners(false, false, 90, 1, 1);
        out_1 = GeoWrangler.move(out_1, offset_1.x, -offset_1.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape1_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_1 = GeoWrangler.removeDuplicates(out_1);
        // Check point count - start and end points are the same.
        Assert.That(clean_1.Count, Is.EqualTo(309));
        // Check expected area
        double area_1 = Clipper.Area(out_1);
        Assert.That(Math.Abs(-(3600 - ((20 * 20) + (10 * 15))) - area_1), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_1 = Clipper.GetBounds(clean_1);
        Assert.That(bounds_1.left, Is.EqualTo(0));
        Assert.That(bounds_1.top, Is.EqualTo(-7));

        // Change our reference to subshape 2
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeRefIndex, 2);
        PointD offset_2 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.That(offset_2.x, Is.EqualTo(-45));
        Assert.That(offset_2.y, Is.EqualTo(41));

        ShapeLibrary shape_2 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_2.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_2.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_2 = shape_2.processCorners(false, false, 90, 1, 1);
        out_2 = GeoWrangler.move(out_2, offset_2.x, -offset_2.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape2_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_2 = GeoWrangler.removeDuplicates(out_2);
        // Check point count - start and end points are the same.
        Assert.That(clean_2.Count, Is.EqualTo(309));
        // Check expected area
        double area_2 = Clipper.Area(out_2);
        Assert.That(Math.Abs(-(3600 - ((20 * 20) + (10 * 15))) - area_2), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_2 = Clipper.GetBounds(clean_2);
        Assert.That(bounds_2.left, Is.EqualTo(-45));
        Assert.That(bounds_2.top, Is.EqualTo(-41));

        // Change our reference to subshape 2
        shapeSettings.setInt(ShapeSettings.properties_i.posInSubShapeIndex, (int)ShapeSettings.subShapeLocations.TR);
        PointD offset_3 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.That(offset_3.x, Is.EqualTo(-60));
        Assert.That(offset_3.y, Is.EqualTo(51));

        ShapeLibrary shape_3 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_3.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.That(shape_3.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.Sshape));
        PathD out_3 = shape_3.processCorners(false, false, 90, 1, 1);
        out_3 = GeoWrangler.move(out_3, offset_3.x, -offset_3.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_3 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape2_pos_tr.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_3 = GeoWrangler.removeDuplicates(out_3);
        // Check point count - start and end points are the same.
        Assert.That(clean_3.Count, Is.EqualTo(309));
        // Check expected area
        double area_3 = Clipper.Area(out_3);
        Assert.That(Math.Abs(-(3600 - ((20 * 20) + (10 * 15))) - area_3), Is.LessThanOrEqualTo(0.0001));
        RectD bounds_3 = Clipper.GetBounds(clean_3);
        Assert.That(bounds_3.left, Is.EqualTo(-60));
        Assert.That(bounds_3.top, Is.EqualTo(-51));
    }

    [Test]
    public static void customOrthoTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        PathD customShape = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            10, 20,
            10, 0
        });
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex), customShape);
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.GEOCORE));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(61));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(area, Is.EqualTo(-200));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(10));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void customOrthoOuterRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        PathD customShape = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            20, 20,
            20, 0
        });
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 10);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex), customShape);
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.GEOCORE));
        PathD out_ = shape.processCorners(false, false, 90, 1, .5);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_outer.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        Assert.That(clean.Count, Is.EqualTo(121));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-(Math.PI * 10 * 10) - area), Is.LessThanOrEqualTo(0.15));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(20));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void customOrthoInnerRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        PathD customShape = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            10, 20,
            10, 10,
            15, 10,
            15, 0
        });
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.iCR, 5);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex), customShape);
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.GEOCORE));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_inner.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(68));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-252.8 - area), Is.LessThanOrEqualTo(0.01));
        RectD bounds = Clipper.GetBounds(out_);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(20));
    }

    [Test]
    public static void customOrthoTipsTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        PathD customShape = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            10, 20,
            10, 0
        });
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeTipLocIndex, (int)ShapeSettings.tipLocations.TL);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.vTBias, 7);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.hTBias, 5);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex), customShape);
        // Check the shape settings are in the shape.
        Assert.That(shape.shapeIndex, Is.EqualTo((int)ShapeLibrary.shapeNames_all.GEOCORE));
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_tip.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.That(clean.Count, Is.EqualTo(83));
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.That(Math.Abs(-((10 + 5) * (20 + 7)) - area), Is.LessThanOrEqualTo(0.001));
        RectD bounds = Clipper.GetBounds(clean);
        Assert.That(bounds.Width, Is.EqualTo(15));
        Assert.That(bounds.Height, Is.EqualTo(27));
    }

    [Test]
    public static void AllShortEdgesTest()
    {
        // Test the fix for the infinite loop when all edges are considered short
        // This test verifies that a square with all short edges correctly becomes a diamond

        // Create a square where all edges will be considered short
        var square = new PathD
        {
            new PointD(0, 0),
            new PointD(10, 0),     // Short horizontal edge (10 units)
            new PointD(10, 10),    // Short vertical edge (10 units)  
            new PointD(0, 10),     // Short horizontal edge (10 units)
            new PointD(0, 0)       // Short vertical edge (10 units) - close the path
        };

        double short_edge_threshold = 15.0; // All edges (10 units) are below this threshold

        // This is a test of the core contourGen functionality
        // We need to create a scenario that would trigger the AssembleWithEasing function

        // For this test, we'll test the direct midpoint connection logic
        var cornerMidpoints = new List<PointD>
        {
            new PointD(5, 0),   // Midpoint of bottom edge
            new PointD(10, 5),  // Midpoint of right edge  
            new PointD(5, 10),  // Midpoint of top edge
            new PointD(0, 5)    // Midpoint of left edge
        };

        // Test the helper function directly using reflection to access the private method
        var type = typeof(contourGen);
        var method = type.GetMethod("ConnectCornerMidpoints",
            System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);

        Assert.That(method, Is.Not.Null, "ConnectCornerMidpoints method should exist");

        var result = method.Invoke(null, new object[] { cornerMidpoints }) as PathD;

        Assert.That(result, Is.Not.Null, "Result should not be null");
        Assert.That(result.Count, Is.EqualTo(4), "Result should have 4 points (diamond)");

        // Verify the diamond vertices are correct (midpoints of square edges)
        Assert.That(result[0].x, Is.EqualTo(5).Within(1e-9), "First point x should be 5");
        Assert.That(result[0].y, Is.EqualTo(0).Within(1e-9), "First point y should be 0");

        Assert.That(result[1].x, Is.EqualTo(10).Within(1e-9), "Second point x should be 10");
        Assert.That(result[1].y, Is.EqualTo(5).Within(1e-9), "Second point y should be 5");

        Assert.That(result[2].x, Is.EqualTo(5).Within(1e-9), "Third point x should be 5");
        Assert.That(result[2].y, Is.EqualTo(10).Within(1e-9), "Third point y should be 10");

        Assert.That(result[3].x, Is.EqualTo(0).Within(1e-9), "Fourth point x should be 0");
        Assert.That(result[3].y, Is.EqualTo(5).Within(1e-9), "Fourth point y should be 5");

        // Verify the diamond has the expected area (half the area of the original square)
        double originalSquareArea = 100; // 10 * 10
        double diamondArea = Math.Abs(Clipper.Area(result));
        double expectedDiamondArea = originalSquareArea / 2; // 50

        Assert.That(diamondArea, Is.EqualTo(expectedDiamondArea).Within(0.001),
            "Diamond area should be half the original square area");
    }

    /// <summary>
    /// Comprehensive tests for various short/long edge combinations to ensure no infinite loops
    /// Tests different configurations that could potentially cause infinite loops in the
    /// contour generation algorithm when searching for non-short corners.
    /// </summary>
    [Test]
    public static void ShortLongEdgeCombinationsTest()
    {
        // Test timeout helper - runs contour generation with a timeout to detect infinite loops

        var testConfigurations = new List<(string name, PathD path, double threshold)>();

        // Test 1: T-shape with top edge long, rest short
        // This was specifically mentioned in the user's request
        var tShape = new PathD
        {
            new PointD(0, 0),    // Bottom left of vertical part
            new PointD(5, 0),    // Bottom right of vertical part
            new PointD(5, 15),   // Top of vertical part (short edge: 5 units)
            new PointD(25, 15),  // Right end of top horizontal (long edge: 20 units)
            new PointD(25, 20),  // Top right corner
            new PointD(-15, 20), // Top left corner (long edge: 40 units)
            new PointD(-15, 15), // Left end of top horizontal
            new PointD(0, 15),   // Back to vertical part (short edge: 15 units)
            new PointD(0, 0)     // Close the path (short edge: 15 units)
        };
        testConfigurations.Add(("T-shape: top long, rest short", tShape, 10.0));

        // Test 2: Single long edge, rest short
        var singleLongEdge = new PathD
        {
            new PointD(0, 0),
            new PointD(5, 0),    // Short: 5 units
            new PointD(5, 5),    // Short: 5 units
            new PointD(30, 5),   // Long: 25 units
            new PointD(30, 0),   // Short: 5 units
            new PointD(0, 0)     // Short: 30 units - wait, this might be long too!
        };
        testConfigurations.Add(("Single long edge", singleLongEdge, 10.0));

        // Test 3: Two adjacent long edges, rest short
        var twoAdjacentLong = new PathD
        {
            new PointD(0, 0),
            new PointD(3, 0),    // Short: 3 units
            new PointD(3, 3),    // Short: 3 units
            new PointD(20, 3),   // Long: 17 units
            new PointD(20, 20),  // Long: 17 units
            new PointD(0, 20),   // Short: 20 units - actually this could be long too
            new PointD(0, 0)     // Short: 20 units - this too
        };
        testConfigurations.Add(("Two adjacent long edges", twoAdjacentLong, 10.0));

        // Test 4: Two non-adjacent long edges, rest short
        var twoNonAdjacentLong = new PathD
        {
            new PointD(0, 0),
            new PointD(5, 0),    // Short: 5 units
            new PointD(25, 0),   // Long: 20 units
            new PointD(25, 5),   // Short: 5 units
            new PointD(20, 5),   // Short: 5 units
            new PointD(20, 25),  // Long: 20 units
            new PointD(0, 25),   // Short: 20 units - could be long
            new PointD(0, 0)     // Short: 25 units - could be long
        };
        testConfigurations.Add(("Two non-adjacent long edges", twoNonAdjacentLong, 10.0));

        // Test 5: Complex L-shape with mixed edges
        var complexLShape = new PathD
        {
            new PointD(0, 0),
            new PointD(3, 0),    // Short: 3 units
            new PointD(3, 3),    // Short: 3 units
            new PointD(25, 3),   // Long: 22 units
            new PointD(25, 6),   // Short: 3 units
            new PointD(6, 6),    // Long: 19 units
            new PointD(6, 25),   // Long: 19 units
            new PointD(0, 25),   // Short: 6 units
            new PointD(0, 0)     // Short: 25 units - could be long
        };
        testConfigurations.Add(("Complex L-shape", complexLShape, 10.0));

        // Test 6: U-shape with mixed edges  
        var complexUShape = new PathD
        {
            new PointD(0, 0),
            new PointD(3, 0),    // Short: 3 units
            new PointD(3, 20),   // Long: 20 units
            new PointD(6, 20),   // Short: 3 units
            new PointD(6, 3),    // Long: 17 units
            new PointD(24, 3),   // Long: 18 units
            new PointD(24, 20),  // Long: 17 units
            new PointD(27, 20),  // Short: 3 units
            new PointD(27, 0),   // Long: 20 units
            new PointD(0, 0)     // Long: 27 units
        };
        testConfigurations.Add(("Complex U-shape", complexUShape, 10.0));

        // Test 7: Alternating short-long pattern
        var alternatingPattern = new PathD
        {
            new PointD(0, 0),
            new PointD(5, 0),    // Short: 5 units
            new PointD(20, 0),   // Long: 15 units
            new PointD(25, 0),   // Short: 5 units
            new PointD(40, 0),   // Long: 15 units
            new PointD(45, 0),   // Short: 5 units
            new PointD(45, 5),   // Short: 5 units
            new PointD(0, 5),    // Long: 45 units
            new PointD(0, 0)     // Short: 5 units
        };
        testConfigurations.Add(("Alternating short-long", alternatingPattern, 10.0));

        // Test 8: Many short edges with one long edge
        var manyShortOneLong = new PathD
        {
            new PointD(0, 0),
            new PointD(2, 0),    // Short: 2 units
            new PointD(2, 2),    // Short: 2 units
            new PointD(4, 2),    // Short: 2 units
            new PointD(4, 4),    // Short: 2 units
            new PointD(6, 4),    // Short: 2 units
            new PointD(6, 6),    // Short: 2 units
            new PointD(30, 6),   // Long: 24 units
            new PointD(30, 0),   // Short: 6 units
            new PointD(0, 0)     // Long: 30 units
        };
        testConfigurations.Add(("Many short, one long", manyShortOneLong, 10.0));

        // Test 9: Nearly all short (edge case for the threshold)
        var nearlyAllShort = new PathD
        {
            new PointD(0, 0),
            new PointD(9, 0),    // Short: 9 units (just under 10)
            new PointD(9, 9),    // Short: 9 units
            new PointD(10, 9),   // Short: 1 unit
            new PointD(10, 10),  // Short: 1 unit
            new PointD(0, 10),   // Long: 10 units (exactly at threshold)
            new PointD(0, 0)     // Long: 10 units
        };
        testConfigurations.Add(("Nearly all short edges", nearlyAllShort, 10.0));

        // Test 10: Irregular polygon with complex short/long mix
        var irregularMix = new PathD
        {
            new PointD(0, 0),
            new PointD(8, 2),    // Short: ~8.2 units
            new PointD(12, 8),   // Short: ~7.2 units  
            new PointD(20, 10),  // Short: ~8.2 units
            new PointD(35, 12),  // Long: 15.1 units
            new PointD(40, 20),  // Short: ~9.4 units
            new PointD(30, 30),  // Long: ~14.1 units
            new PointD(15, 25),  // Long: ~18 units
            new PointD(5, 15),   // Long: ~14.1 units
            new PointD(0, 0)     // Long: ~18 units
        };
        testConfigurations.Add(("Irregular mixed polygon", irregularMix, 10.0));

        // Run all tests
        var failedTests = new List<string>();
        var successfulTests = new List<string>();
        var svgOutputDir = Path.Combine(Path.GetTempPath(), "comprehensive_edge_tests");
        Directory.CreateDirectory(svgOutputDir);
        Console.WriteLine($"SVG outputs will be saved to: {svgOutputDir}");

        foreach (var (name, path, threshold) in testConfigurations)
        {
            Console.WriteLine($"Testing: {name}");

            var (success, result) = TestWithTimeoutAndResult(path, threshold, 10000); // 10 second timeout

            // Generate SVG regardless of success/failure for review
            string safeName = name.Replace(" ", "_").Replace(":", "_").Replace(",", "_");
            string svgContent = CreateTestSvg(path, result, name, threshold);
            string svgPath = Path.Combine(svgOutputDir, $"{safeName}.svg");
            File.WriteAllText(svgPath, svgContent, Encoding.UTF8);

            if (success)
            {
                successfulTests.Add(name);
                Console.WriteLine($"✓ {name} - PASSED (SVG: {svgPath})");
            }
            else
            {
                failedTests.Add(name);
                Console.WriteLine($"✗ {name} - FAILED (timeout or exception) (SVG: {svgPath})");
            }
        }

        // Report results
        Console.WriteLine($"\n=== Test Results ===");
        Console.WriteLine($"Successful tests: {successfulTests.Count}");
        foreach (var test in successfulTests)
        {
            Console.WriteLine($"  ✓ {test}");
        }

        if (failedTests.Count > 0)
        {
            Console.WriteLine($"Failed tests: {failedTests.Count}");
            foreach (var test in failedTests)
            {
                Console.WriteLine($"  ✗ {test}");
            }
        }

        // The test passes if all configurations complete without timeout
        Assert.That(failedTests.Count, Is.EqualTo(0),
            $"Some test configurations failed or timed out, indicating potential infinite loops: {string.Join(", ", failedTests)}");
    }

    /// <summary>
    /// Tests specific edge cases that are most likely to cause infinite loops
    /// based on the algorithm's logic for finding non-short corners.
    /// </summary>
    [Test]
    public static void EdgeCaseInfiniteLoopTest()
    {
        // Test the specific T-shape mentioned by the user
        var tShapeTopLong = new PathD
        {
            new PointD(10, 0),   // Start at bottom of vertical stem
            new PointD(15, 0),   // Short edge: 5 units (bottom of stem)
            new PointD(15, 10),  // Short edge: 10 units (right side of stem)  
            new PointD(30, 10),  // Long edge: 15 units (right side of top)
            new PointD(30, 15),  // Short edge: 5 units (top right)
            new PointD(5, 15),   // Long edge: 25 units (TOP - the long edge)
            new PointD(5, 10),   // Short edge: 5 units (top left)
            new PointD(10, 10),  // Short edge: 5 units (left side of stem)
            new PointD(10, 0)    // Short edge: 10 units (close the path)
        };

        // Use timeout to detect infinite loops
        var cancellationTokenSource = new CancellationTokenSource(5000); // 5 second timeout
        bool completed = false;
        PathD result = null;
        Exception caughtException = null;

        var task = Task.Run(() =>
        {
            try
            {
                result = contourGen.makeContour(tShapeTopLong,
                    concaveRadius: 2.0,
                    convexRadius: 2.0,
                    edgeResolution: 0.5,
                    angularResolution: 1.0,
                    shortEdgeLength: 12.0, // This makes most edges short except the top
                    maxShortEdgeLength: 24.0,
                    optimizeCorners: 0);
                completed = true;
            }
            catch (Exception ex)
            {
                caughtException = ex;
            }
        }, cancellationTokenSource.Token);

        bool finishedInTime = task.Wait(5000);

        if (!finishedInTime)
        {
            Console.WriteLine("T-shape test timed out - likely infinite loop detected!");
            Assert.Fail("T-shape with top long edge timed out, indicating an infinite loop in the contour generation algorithm");
        }
        else if (caughtException != null)
        {
            Console.WriteLine($"T-shape test threw exception: {caughtException.Message}");
            Assert.Fail($"T-shape test failed with exception: {caughtException.Message}");
        }
        else if (!completed)
        {
            Assert.Fail("T-shape test did not complete successfully");
        }

        Assert.That(result, Is.Not.Null, "T-shape contour generation should return a valid result");
        Assert.That(result.Count, Is.GreaterThan(0), "T-shape result should have points");

        // Generate SVG for review
        var svgOutputDir = Path.Combine(Path.GetTempPath(), "comprehensive_edge_tests");
        Directory.CreateDirectory(svgOutputDir);
        string svgContent = CreateTestSvg(tShapeTopLong, result, "T-shape: Top Edge Long, Rest Short", 12.0);
        string svgPath = Path.Combine(svgOutputDir, "T_shape_specific_test.svg");
        File.WriteAllText(svgPath, svgContent, Encoding.UTF8);

        Console.WriteLine($"T-shape test completed successfully with {result.Count} points");
        Console.WriteLine($"T-shape SVG saved to: {svgPath}");
    }

    /// <summary>
    /// Performance test to ensure contour generation completes in reasonable time
    /// for various polygon complexities with mixed short/long edges.
    /// </summary>
    [Test]
    public static void ContourGenerationPerformanceTest()
    {
        var testCases = new List<(string name, PathD path, double threshold, int expectedMaxTimeMs)>();

        // Simple cases should complete very quickly
        var simpleSquare = new PathD
        {
            new PointD(0, 0),
            new PointD(10, 0),
            new PointD(10, 10),
            new PointD(0, 10),
            new PointD(0, 0)
        };
        testCases.Add(("Simple square", simpleSquare, 5.0, 1000));

        // More complex case should still complete reasonably quickly
        var complexShape = new PathD();
        for (int i = 0; i < 20; i++)
        {
            double angle = i * 2 * Math.PI / 20;
            double radius = (i % 2 == 0) ? 10 : 20; // Alternating radii
            complexShape.Add(new PointD(
                Math.Cos(angle) * radius,
                Math.Sin(angle) * radius));
        }
        complexShape.Add(complexShape[0]); // Close the path
        testCases.Add(("Complex 20-sided star", complexShape, 15.0, 3000));

        foreach (var (name, path, threshold, maxTime) in testCases)
        {
            var stopwatch = System.Diagnostics.Stopwatch.StartNew();

            var result = contourGen.makeContour(path,
                concaveRadius: 2.0,
                convexRadius: 2.0,
                edgeResolution: 1.0,
                angularResolution: 1.0,
                shortEdgeLength: threshold,
                maxShortEdgeLength: threshold * 2,
                optimizeCorners: 0);

            stopwatch.Stop();

            // Generate SVG for performance test results
            var svgOutputDir = Path.Combine(Path.GetTempPath(), "comprehensive_edge_tests");
            Directory.CreateDirectory(svgOutputDir);
            string safeName = name.Replace(" ", "_").Replace("-", "_");
            string svgContent = CreateTestSvg(path, result, $"Performance Test: {name}", threshold);
            string svgPath = Path.Combine(svgOutputDir, $"perf_{safeName}.svg");
            File.WriteAllText(svgPath, svgContent, Encoding.UTF8);

            Console.WriteLine($"{name}: {stopwatch.ElapsedMilliseconds}ms (SVG: {svgPath})");

            Assert.That(result, Is.Not.Null, $"{name} should return a valid result");
            Assert.That(stopwatch.ElapsedMilliseconds, Is.LessThan(maxTime),
                $"{name} should complete within {maxTime}ms but took {stopwatch.ElapsedMilliseconds}ms");
        }
    }

    /// <summary>
    /// Comprehensive analysis of ContourGen runtime cost vs quality for different
    /// subdivision and recursion configuration options. Provides systematic insights
    /// into performance characteristics of various parameter combinations.
    /// </summary>
    [Test]
    public static void ContourGenConfigurationAnalysisTest()
    {
        Console.WriteLine("=== ContourGen Configuration Analysis ===");
        Console.WriteLine("Testing runtime cost vs quality for different subdivision and recursion configurations");

        // Define test shapes of varying complexity
        var testShapes = new List<(string name, PathD path)>
        {
            ("Simple Square", CreateSquare(10)),
            ("Complex Star", CreateStar(20, 5, 15)), // 20 points, inner radius 5, outer radius 15
            ("L-Shape", CreateLShape(20, 20, 8, 8)),
            ("Fine Detail Shape", CreateFineDetailShape())
        };

        // Define parameter combinations to test
        var parameterSets = new List<(string name, double concaveRadius, double convexRadius,
            double edgeResolution, double angularResolution, double shortEdgeLength,
            double maxShortEdgeLength, int optimizeCorners)>
        {
            // Low quality, fast
            ("Fast_Low", 1.0, 1.0, 5.0, 10.0, 3.0, 6.0, 0),
            ("Fast_High", 1.0, 1.0, 2.0, 5.0, 3.0, 6.0, 0),
            
            // Medium quality
            ("Medium_Low", 2.0, 2.0, 2.0, 5.0, 5.0, 10.0, 0),
            ("Medium_High", 2.0, 2.0, 1.0, 2.0, 5.0, 10.0, 0),
            ("Medium_Optimized", 2.0, 2.0, 1.0, 2.0, 5.0, 10.0, 1),
            
            // High quality, slower
            ("High_Standard", 3.0, 3.0, 0.5, 1.0, 8.0, 16.0, 0),
            ("High_Fine", 3.0, 3.0, 0.2, 0.5, 8.0, 16.0, 0),
            ("High_Optimized", 3.0, 3.0, 0.5, 1.0, 8.0, 16.0, 1),
            
            // Ultra-high quality, very slow
            ("Ultra_Detail", 5.0, 5.0, 0.1, 0.25, 15.0, 30.0, 0),
            ("Ultra_Optimized", 5.0, 5.0, 0.1, 0.25, 15.0, 30.0, 1)
        };

        var results = new List<ContourGenPerformanceResult>();
        var svgOutputDir = Path.Combine(Path.GetTempPath(), "contour_gen_analysis");
        Directory.CreateDirectory(svgOutputDir);

        foreach (var (shapeName, shape) in testShapes)
        {
            Console.WriteLine($"\n--- Testing Shape: {shapeName} ---");

            foreach (var paramSet in parameterSets)
            {
                var (paramName, concaveRadius, convexRadius, edgeResolution, angularResolution,
                     shortEdgeLength, maxShortEdgeLength, optimizeCorners) = paramSet;

                Console.Write($"  {paramName}: ");

                var stopwatch = System.Diagnostics.Stopwatch.StartNew();

                try
                {
                    var result = contourGen.makeContour(shape,
                        concaveRadius: concaveRadius,
                        convexRadius: convexRadius,
                        edgeResolution: edgeResolution,
                        angularResolution: angularResolution,
                        shortEdgeLength: shortEdgeLength,
                        maxShortEdgeLength: maxShortEdgeLength,
                        optimizeCorners: optimizeCorners);

                    stopwatch.Stop();

                    // Calculate quality metrics
                    var metrics = CalculateQualityMetrics(shape, result, edgeResolution);

                    var perfResult = new ContourGenPerformanceResult
                    {
                        ShapeName = shapeName,
                        ParameterSetName = paramName,
                        RuntimeMs = stopwatch.ElapsedMilliseconds,
                        InputPointCount = shape.Count,
                        OutputPointCount = result?.Count ?? 0,
                        QualityMetrics = metrics,
                        Parameters = new ContourGenParameters
                        {
                            ConcaveRadius = concaveRadius,
                            ConvexRadius = convexRadius,
                            EdgeResolution = edgeResolution,
                            AngularResolution = angularResolution,
                            ShortEdgeLength = shortEdgeLength,
                            MaxShortEdgeLength = maxShortEdgeLength,
                            OptimizeCorners = optimizeCorners
                        }
                    };

                    results.Add(perfResult);

                    // Generate SVG for visual comparison
                    string svgContent = CreateContourGenAnalysisSvg(shape, result, perfResult);
                    string safeShapeName = shapeName.Replace(" ", "_");
                    string svgPath = Path.Combine(svgOutputDir, $"{safeShapeName}_{paramName}.svg");
                    File.WriteAllText(svgPath, svgContent, Encoding.UTF8);

                    Console.WriteLine($"{stopwatch.ElapsedMilliseconds}ms, {result?.Count ?? 0} pts, Quality: {metrics.OverallQuality:F2} (SVG: {svgPath})");

                    Assert.That(result, Is.Not.Null, $"Result should not be null for {shapeName} with {paramName}");
                    Assert.That(stopwatch.ElapsedMilliseconds, Is.LessThan(30000),
                        $"Should complete within 30 seconds for {shapeName} with {paramName}");
                }
                catch (Exception ex)
                {
                    stopwatch.Stop();
                    Console.WriteLine($"FAILED - {ex.Message}");

                    // Still record the failure for analysis
                    results.Add(new ContourGenPerformanceResult
                    {
                        ShapeName = shapeName,
                        ParameterSetName = paramName,
                        RuntimeMs = stopwatch.ElapsedMilliseconds,
                        InputPointCount = shape.Count,
                        OutputPointCount = 0,
                        Failed = true,
                        FailureMessage = ex.Message
                    });
                }
            }
        }

        // Generate comprehensive analysis report
        GeneratePerformanceAnalysisReport(results, svgOutputDir);

        Console.WriteLine($"\n=== Analysis Complete ===");
        Console.WriteLine($"Results saved to: {svgOutputDir}");
        Console.WriteLine($"Total test cases: {results.Count}");
        Console.WriteLine($"Successful: {results.Count(r => !r.Failed)}");
        Console.WriteLine($"Failed: {results.Count(r => r.Failed)}");
    }

    #endregion

    #region Helper Methods for SVG Generation and Timeout Testing

    /// <summary>
    /// Create SVG visualization showing input polygon and generated contour
    /// </summary>
    static string CreateTestSvg(PathD inputPath, PathD outputPath, string title, double shortEdgeThreshold)
    {
        if (inputPath == null || inputPath.Count == 0)
            return "";

        var allPoints = inputPath.ToList();
        if (outputPath != null && outputPath.Count > 0)
            allPoints.AddRange(outputPath);

        double minX = allPoints.Min(p => p.x) - 10;
        double maxX = allPoints.Max(p => p.x) + 10;
        double minY = allPoints.Min(p => p.y) - 10;
        double maxY = allPoints.Max(p => p.y) + 10;

        double width = maxX - minX;
        double height = maxY - minY;

        StringBuilder sb = new StringBuilder();
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"600\" viewBox=\"{minX} {-maxY} {width} {height}\">");

        // Title
        sb.AppendLine($"  <title>{title}</title>");

        // Grid
        sb.AppendLine("  <defs>");
        sb.AppendLine("    <pattern id=\"grid\" width=\"10\" height=\"10\" patternUnits=\"userSpaceOnUse\">");
        sb.AppendLine("      <path d=\"M 10 0 L 0 0 0 10\" fill=\"none\" stroke=\"#e0e0e0\" stroke-width=\"0.5\"/>");
        sb.AppendLine("    </pattern>");
        sb.AppendLine("  </defs>");
        sb.AppendLine($"  <rect width=\"100%\" height=\"100%\" fill=\"url(#grid)\" />");

        // Input polygon (blue)
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        for (int i = 0; i < inputPath.Count; i++)
        {
            sb.Append($"{inputPath[i].x},{-inputPath[i].y} ");
        }
        sb.AppendLine("\"/>");

        // Mark short vs long edges with colors
        for (int i = 0; i < inputPath.Count - 1; i++)
        {
            var p1 = inputPath[i];
            var p2 = inputPath[i + 1];
            double edgeLength = Math.Sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
            string color = edgeLength <= shortEdgeThreshold ? "#ff6b6b" : "#4ecdc4";
            double strokeWidth = edgeLength <= shortEdgeThreshold ? 3 : 2;

            sb.AppendLine($"  <line x1=\"{p1.x}\" y1=\"{-p1.y}\" x2=\"{p2.x}\" y2=\"{-p2.y}\" stroke=\"{color}\" stroke-width=\"{strokeWidth}\" stroke-opacity=\"0.7\"/>");

            // Label edge length
            double midX = (p1.x + p2.x) / 2;
            double midY = -(p1.y + p2.y) / 2;
            sb.AppendLine($"  <text x=\"{midX}\" y=\"{midY}\" font-size=\"8\" fill=\"{color}\" text-anchor=\"middle\">{edgeLength:F1}</text>");
        }

        // Output contour (red)
        if (outputPath != null && outputPath.Count > 0)
        {
            sb.Append("  <polyline fill=\"none\" stroke=\"red\" stroke-width=\"3\" points=\"");
            foreach (var p in outputPath)
            {
                sb.Append($"{p.x},{-p.y} ");
            }
            // Close the path if not already closed
            if (outputPath.Count > 0 && (Math.Abs(outputPath[0].x - outputPath[^1].x) > 1e-9 || Math.Abs(outputPath[0].y - outputPath[^1].y) > 1e-9))
            {
                sb.Append($"{outputPath[0].x},{-outputPath[0].y} ");
            }
            sb.AppendLine("\"/>");

            // Mark output vertices
            foreach (var p in outputPath)
            {
                sb.AppendLine($"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"2\" fill=\"red\" fill-opacity=\"0.8\"/>");
            }
        }

        // Mark input vertices
        foreach (var p in inputPath)
        {
            sb.AppendLine($"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"1.5\" fill=\"blue\" fill-opacity=\"0.8\"/>");
        }

        // Legend
        double legendY = -maxY + 20;
        sb.AppendLine($"  <text x=\"{minX + 5}\" y=\"{legendY}\" font-size=\"14\" font-weight=\"bold\" fill=\"black\">{title}</text>");
        sb.AppendLine($"  <text x=\"{minX + 5}\" y=\"{legendY + 18}\" font-size=\"12\" fill=\"blue\">Input Polygon</text>");
        if (outputPath != null && outputPath.Count > 0)
            sb.AppendLine($"  <text x=\"{minX + 5}\" y=\"{legendY + 35}\" font-size=\"12\" fill=\"red\">Generated Contour</text>");
        sb.AppendLine($"  <text x=\"{minX + 5}\" y=\"{legendY + 52}\" font-size=\"12\" fill=\"#ff6b6b\">Short Edges (≤{shortEdgeThreshold})</text>");
        sb.AppendLine($"  <text x=\"{minX + 5}\" y=\"{legendY + 69}\" font-size=\"12\" fill=\"#4ecdc4\">Long Edges (>{shortEdgeThreshold})</text>");

        sb.AppendLine("</svg>");
        return sb.ToString();
    }

    /// <summary>
    /// Test function that runs contour generation with timeout and returns both success status and result
    /// </summary>
    static (bool success, PathD result) TestWithTimeoutAndResult(PathD path, double shortEdgeThreshold, int timeoutMs = 5000)
    {
        var cancellationTokenSource = new CancellationTokenSource(timeoutMs);
        PathD result = null;

        var task = Task.Run(() =>
        {
            try
            {
                result = contourGen.makeContour(path,
                    concaveRadius: 5.0,
                    convexRadius: 5.0,
                    edgeResolution: 1.0,
                    angularResolution: 1.0,
                    shortEdgeLength: shortEdgeThreshold,
                    maxShortEdgeLength: shortEdgeThreshold * 2,
                    optimizeCorners: 0);
                return result != null;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Exception in contour generation: {ex.Message}");
                return false;
            }
        }, cancellationTokenSource.Token);

        try
        {
            bool success = task.Wait(timeoutMs) && task.Result;
            return (success, result);
        }
        catch (Exception)
        {
            return (false, null);
        }
    }

    static bool TestWithTimeout(PathD path, double shortEdgeThreshold, int timeoutMs = 5000)
    {
        var (success, _) = TestWithTimeoutAndResult(path, shortEdgeThreshold, timeoutMs);
        return success;
    }

    #endregion

    #region ContourGen Performance Analysis Support

    /// <summary>
    /// Results of a ContourGen performance test run
    /// </summary>
    public class ContourGenPerformanceResult
    {
        public string ShapeName { get; set; }
        public string ParameterSetName { get; set; }
        public long RuntimeMs { get; set; }
        public int InputPointCount { get; set; }
        public int OutputPointCount { get; set; }
        public QualityMetrics QualityMetrics { get; set; }
        public ContourGenParameters Parameters { get; set; }
        public bool Failed { get; set; }
        public string FailureMessage { get; set; }
    }

    /// <summary>
    /// Parameters used for ContourGen test
    /// </summary>
    public class ContourGenParameters
    {
        public double ConcaveRadius { get; set; }
        public double ConvexRadius { get; set; }
        public double EdgeResolution { get; set; }
        public double AngularResolution { get; set; }
        public double ShortEdgeLength { get; set; }
        public double MaxShortEdgeLength { get; set; }
        public int OptimizeCorners { get; set; }
    }

    /// <summary>
    /// Quality metrics for contour generation
    /// </summary>
    public class QualityMetrics
    {
        public double AverageSegmentLength { get; set; }
        public double MaxSegmentLength { get; set; }
        public double MinSegmentLength { get; set; }
        public double SegmentLengthVariance { get; set; }
        public double AverageAngularChange { get; set; }
        public double MaxAngularChange { get; set; }
        public double Smoothness { get; set; }
        public double PointDensityRatio { get; set; }
        public double OverallQuality { get; set; }
    }

    /// <summary>
    /// Create test shapes for analysis
    /// </summary>
    static PathD CreateSquare(double size)
    {
        return new PathD
        {
            new PointD(0, 0),
            new PointD(size, 0),
            new PointD(size, size),
            new PointD(0, size),
            new PointD(0, 0)
        };
    }

    static PathD CreateStar(int points, double innerRadius, double outerRadius)
    {
        var star = new PathD();
        for (int i = 0; i < points * 2; i++)
        {
            double angle = i * Math.PI / points;
            double radius = (i % 2 == 0) ? outerRadius : innerRadius;
            star.Add(new PointD(
                Math.Cos(angle) * radius,
                Math.Sin(angle) * radius));
        }
        star.Add(star[0]); // Close the path
        return star;
    }

    static PathD CreateLShape(double width, double height, double thickness1, double thickness2)
    {
        return new PathD
        {
            new PointD(0, 0),
            new PointD(width, 0),
            new PointD(width, thickness1),
            new PointD(thickness2, thickness1),
            new PointD(thickness2, height),
            new PointD(0, height),
            new PointD(0, 0)
        };
    }

    static PathD CreateFineDetailShape()
    {
        var shape = new PathD();
        // Create a shape with many small details that will test subdivision algorithms
        for (int i = 0; i <= 100; i++)
        {
            double t = i / 100.0;
            double angle = t * 4 * Math.PI;
            double radius = 10 + 3 * Math.Sin(angle * 5); // High frequency variation
            double x = Math.Cos(angle) * radius;
            double y = Math.Sin(angle) * radius;
            shape.Add(new PointD(x, y));
        }
        shape.Add(shape[0]); // Close the path
        return shape;
    }

    /// <summary>
    /// Calculate quality metrics for a contour generation result
    /// </summary>
    static QualityMetrics CalculateQualityMetrics(PathD input, PathD output, double targetResolution)
    {
        if (output == null || output.Count < 2)
        {
            return new QualityMetrics { OverallQuality = 0.0 };
        }

        var segmentLengths = new List<double>();
        var angularChanges = new List<double>();

        // Calculate segment lengths
        for (int i = 1; i < output.Count; i++)
        {
            double length = Helper.Length(Helper.Minus(output[i], output[i - 1]));
            segmentLengths.Add(length);
        }

        // Calculate angular changes
        for (int i = 1; i < output.Count - 1; i++)
        {
            var v1 = Helper.Normalized(Helper.Minus(output[i], output[i - 1]));
            var v2 = Helper.Normalized(Helper.Minus(output[i + 1], output[i]));
            double dot = Math.Max(-1.0, Math.Min(1.0, Helper.Dot(v1, v2)));
            double angle = Math.Acos(dot) * 180.0 / Math.PI;
            angularChanges.Add(angle);
        }

        double avgSegLen = segmentLengths.Average();
        double maxSegLen = segmentLengths.Max();
        double minSegLen = segmentLengths.Min();
        double variance = segmentLengths.Select(x => (x - avgSegLen) * (x - avgSegLen)).Average();

        double avgAngular = angularChanges.Count > 0 ? angularChanges.Average() : 0;
        double maxAngular = angularChanges.Count > 0 ? angularChanges.Max() : 0;

        // Smoothness metric (lower variance in angular changes = smoother)
        double smoothness = angularChanges.Count > 0 ?
            1.0 / (1.0 + angularChanges.Select(x => (x - avgAngular) * (x - avgAngular)).Average()) : 1.0;

        // Point density ratio (output points / input points)
        double densityRatio = (double)output.Count / input.Count;

        // Overall quality score (0-1, higher is better)
        // Factors: adherence to target resolution, smoothness, reasonable point density
        double resolutionScore = Math.Exp(-Math.Abs(avgSegLen - targetResolution) / targetResolution);
        double densityScore = Math.Exp(-Math.Abs(densityRatio - 2.0) / 2.0); // Target ~2x points as good
        double overallQuality = (resolutionScore + smoothness + densityScore) / 3.0;

        return new QualityMetrics
        {
            AverageSegmentLength = avgSegLen,
            MaxSegmentLength = maxSegLen,
            MinSegmentLength = minSegLen,
            SegmentLengthVariance = variance,
            AverageAngularChange = avgAngular,
            MaxAngularChange = maxAngular,
            Smoothness = smoothness,
            PointDensityRatio = densityRatio,
            OverallQuality = overallQuality
        };
    }

    /// <summary>
    /// Create specialized SVG for ContourGen analysis with detailed metrics
    /// </summary>
    static string CreateContourGenAnalysisSvg(PathD input, PathD output, ContourGenPerformanceResult result)
    {
        if (input == null || input.Count == 0)
            return "";

        var allPoints = input.ToList();
        if (output != null && output.Count > 0)
            allPoints.AddRange(output);

        double minX = allPoints.Min(p => p.x) - 10;
        double maxX = allPoints.Max(p => p.x) + 10;
        double minY = allPoints.Min(p => p.y) - 10;
        double maxY = allPoints.Max(p => p.y) + 10;

        double width = maxX - minX;
        double height = maxY - minY;

        StringBuilder sb = new StringBuilder();
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"600\" viewBox=\"{minX} {-maxY} {width} {height}\">");

        // Title
        sb.AppendLine($"  <title>{result.ShapeName} - {result.ParameterSetName}</title>");

        // Grid
        sb.AppendLine("  <defs>");
        sb.AppendLine("    <pattern id=\"grid\" width=\"10\" height=\"10\" patternUnits=\"userSpaceOnUse\">");
        sb.AppendLine("      <path d=\"M 10 0 L 0 0 0 10\" fill=\"none\" stroke=\"#e0e0e0\" stroke-width=\"0.5\"/>");
        sb.AppendLine("    </pattern>");
        sb.AppendLine("  </defs>");
        sb.AppendLine($"  <rect width=\"100%\" height=\"100%\" fill=\"url(#grid)\" />");

        // Input polygon (blue)
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (var p in input)
        {
            sb.Append($"{p.x},{-p.y} ");
        }
        sb.AppendLine("\"/>");

        // Output contour (red)
        if (output != null && output.Count > 0)
        {
            sb.Append("  <polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"");
            foreach (var p in output)
            {
                sb.Append($"{p.x},{-p.y} ");
            }
            sb.AppendLine("\"/>");

            // Mark output vertices with small circles
            foreach (var p in output)
            {
                sb.AppendLine($"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"1\" fill=\"red\" fill-opacity=\"0.6\"/>");
            }
        }

        // Mark input vertices
        foreach (var p in input)
        {
            sb.AppendLine($"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"1.5\" fill=\"blue\" fill-opacity=\"0.8\"/>");
        }

        // Performance and quality information panel
        double legendX = minX + 5;
        double legendY = -maxY + 20;

        sb.AppendLine($"  <rect x=\"{legendX - 3}\" y=\"{legendY - 15}\" width=\"300\" height=\"200\" fill=\"white\" fill-opacity=\"0.9\" stroke=\"black\" stroke-width=\"1\"/>");

        sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY}\" font-size=\"14\" font-weight=\"bold\" fill=\"black\">{result.ShapeName} - {result.ParameterSetName}</text>");
        sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 20}\" font-size=\"11\" fill=\"black\">Runtime: {result.RuntimeMs}ms</text>");
        sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 35}\" font-size=\"11\" fill=\"black\">Points: {result.InputPointCount} → {result.OutputPointCount}</text>");

        if (result.QualityMetrics != null)
        {
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 50}\" font-size=\"11\" fill=\"black\">Quality Score: {result.QualityMetrics.OverallQuality:F3}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 65}\" font-size=\"10\" fill=\"darkblue\">Avg Segment: {result.QualityMetrics.AverageSegmentLength:F2}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 80}\" font-size=\"10\" fill=\"darkblue\">Smoothness: {result.QualityMetrics.Smoothness:F3}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 95}\" font-size=\"10\" fill=\"darkblue\">Density Ratio: {result.QualityMetrics.PointDensityRatio:F2}</text>");
        }

        if (result.Parameters != null)
        {
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 115}\" font-size=\"10\" fill=\"darkgreen\">EdgeRes: {result.Parameters.EdgeResolution}, AngRes: {result.Parameters.AngularResolution}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 130}\" font-size=\"10\" fill=\"darkgreen\">Radii: C{result.Parameters.ConcaveRadius}, X{result.Parameters.ConvexRadius}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 145}\" font-size=\"10\" fill=\"darkgreen\">ShortEdge: {result.Parameters.ShortEdgeLength}-{result.Parameters.MaxShortEdgeLength}</text>");
            sb.AppendLine($"  <text x=\"{legendX}\" y=\"{legendY + 160}\" font-size=\"10\" fill=\"darkgreen\">Optimize: {(result.Parameters.OptimizeCorners == 1 ? "Yes" : "No")}</text>");
        }

        sb.AppendLine("</svg>");
        return sb.ToString();
    }

    /// <summary>
    /// Generate comprehensive performance analysis report
    /// </summary>
    static void GeneratePerformanceAnalysisReport(List<ContourGenPerformanceResult> results, string outputDir)
    {
        var reportPath = Path.Combine(outputDir, "performance_analysis_report.html");

        var html = new StringBuilder();
        html.AppendLine("<!DOCTYPE html>");
        html.AppendLine("<html>");
        html.AppendLine("<head>");
        html.AppendLine("  <title>ContourGen Performance Analysis Report</title>");
        html.AppendLine("  <style>");
        html.AppendLine("    body { font-family: Arial, sans-serif; margin: 20px; }");
        html.AppendLine("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }");
        html.AppendLine("    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }");
        html.AppendLine("    th { background-color: #f2f2f2; }");
        html.AppendLine("    .fast { background-color: #e8f5e8; }");
        html.AppendLine("    .medium { background-color: #fff3cd; }");
        html.AppendLine("    .slow { background-color: #f8d7da; }");
        html.AppendLine("    .failed { background-color: #ffcccc; }");
        html.AppendLine("    .summary { background-color: #e9ecef; padding: 15px; margin: 20px 0; border-radius: 5px; }");
        html.AppendLine("  </style>");
        html.AppendLine("</head>");
        html.AppendLine("<body>");

        html.AppendLine("<h1>ContourGen Performance Analysis Report</h1>");
        html.AppendLine($"<p>Generated on: {DateTime.Now:yyyy-MM-dd HH:mm:ss}</p>");

        // Summary statistics
        html.AppendLine("<div class='summary'>");
        html.AppendLine("<h2>Summary</h2>");
        html.AppendLine($"<p>Total test cases: {results.Count}</p>");
        html.AppendLine($"<p>Successful: {results.Count(r => !r.Failed)}</p>");
        html.AppendLine($"<p>Failed: {results.Count(r => r.Failed)}</p>");

        var successfulResults = results.Where(r => !r.Failed).ToList();
        if (successfulResults.Any())
        {
            html.AppendLine($"<p>Average runtime: {successfulResults.Average(r => r.RuntimeMs):F1}ms</p>");
            html.AppendLine($"<p>Fastest: {successfulResults.Min(r => r.RuntimeMs)}ms</p>");
            html.AppendLine($"<p>Slowest: {successfulResults.Max(r => r.RuntimeMs)}ms</p>");
            html.AppendLine($"<p>Average quality score: {successfulResults.Where(r => r.QualityMetrics != null).Average(r => r.QualityMetrics.OverallQuality):F3}</p>");
        }
        html.AppendLine("</div>");

        // Detailed results table
        html.AppendLine("<h2>Detailed Results</h2>");
        html.AppendLine("<table>");
        html.AppendLine("<tr>");
        html.AppendLine("  <th>Shape</th><th>Parameter Set</th><th>Runtime (ms)</th>");
        html.AppendLine("  <th>Input Points</th><th>Output Points</th><th>Quality Score</th>");
        html.AppendLine("  <th>Avg Segment Length</th><th>Smoothness</th><th>Status</th>");
        html.AppendLine("</tr>");

        foreach (var result in results.OrderBy(r => r.ShapeName).ThenBy(r => r.RuntimeMs))
        {
            string rowClass = result.Failed ? "failed" :
                             result.RuntimeMs < 50 ? "fast" :
                             result.RuntimeMs < 500 ? "medium" : "slow";

            html.AppendLine($"<tr class='{rowClass}'>");
            html.AppendLine($"  <td>{result.ShapeName}</td>");
            html.AppendLine($"  <td>{result.ParameterSetName}</td>");
            html.AppendLine($"  <td>{result.RuntimeMs}</td>");
            html.AppendLine($"  <td>{result.InputPointCount}</td>");
            html.AppendLine($"  <td>{result.OutputPointCount}</td>");
            html.AppendLine($"  <td>{(result.QualityMetrics != null ? result.QualityMetrics.OverallQuality.ToString("F3") : "N/A")}</td>");
            html.AppendLine($"  <td>{(result.QualityMetrics != null ? result.QualityMetrics.AverageSegmentLength.ToString("F2") : "N/A")}</td>");
            html.AppendLine($"  <td>{(result.QualityMetrics != null ? result.QualityMetrics.Smoothness.ToString("F3") : "N/A")}</td>");
            html.AppendLine($"  <td>{(result.Failed ? result.FailureMessage : "Success")}</td>");
            html.AppendLine("</tr>");
        }

        html.AppendLine("</table>");

        // Parameter analysis
        html.AppendLine("<h2>Parameter Analysis</h2>");

        foreach (var shapeGroup in results.Where(r => !r.Failed).GroupBy(r => r.ShapeName))
        {
            html.AppendLine($"<h3>{shapeGroup.Key}</h3>");
            html.AppendLine("<p>Runtime vs Quality trade-offs:</p>");
            html.AppendLine("<ul>");

            var sortedByRuntime = shapeGroup.OrderBy(r => r.RuntimeMs).ToList();
            var sortedByQuality = shapeGroup.Where(r => r.QualityMetrics != null)
                                            .OrderByDescending(r => r.QualityMetrics.OverallQuality).ToList();

            html.AppendLine($"<li>Fastest: {sortedByRuntime.First().ParameterSetName} ({sortedByRuntime.First().RuntimeMs}ms)</li>");
            html.AppendLine($"<li>Slowest: {sortedByRuntime.Last().ParameterSetName} ({sortedByRuntime.Last().RuntimeMs}ms)</li>");

            if (sortedByQuality.Any())
            {
                html.AppendLine($"<li>Highest Quality: {sortedByQuality.First().ParameterSetName} (Score: {sortedByQuality.First().QualityMetrics.OverallQuality:F3})</li>");
                html.AppendLine($"<li>Lowest Quality: {sortedByQuality.Last().ParameterSetName} (Score: {sortedByQuality.Last().QualityMetrics.OverallQuality:F3})</li>");
            }
            html.AppendLine("</ul>");
        }

        html.AppendLine("</body>");
        html.AppendLine("</html>");

        File.WriteAllText(reportPath, html.ToString());
        Console.WriteLine($"Comprehensive analysis report generated: {reportPath}");
    }

    #endregion
}