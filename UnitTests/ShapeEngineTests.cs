using Clipper2Lib;
using geoWrangler;
using shapeEngine;

namespace UnitTests;

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
        Assert.AreEqual(shapeTable.Length, shapes.Count);
        Assert.AreEqual(shapeTable2.Length, shapes2.Count);
        Assert.AreEqual(0, shapes.IndexOf(availableShapes_all[0]));
        Assert.AreEqual(1, shapes.IndexOf(availableShapes_all[1]));
        Assert.AreEqual(2, shapes.IndexOf(availableShapes_all[2]));
        Assert.AreEqual(3, shapes.IndexOf(availableShapes_all[3]));
        Assert.AreEqual(4, shapes.IndexOf(availableShapes_all[4]));
        Assert.AreEqual(5, shapes.IndexOf(availableShapes_all[5]));
        Assert.AreEqual(6, shapes.IndexOf(availableShapes_all[6]));
        Assert.AreEqual(7, shapes.IndexOf(availableShapes_all[7]));
        Assert.AreEqual(8, shapes.IndexOf(availableShapes_all[8]));
        Assert.AreEqual(-1, shapes.IndexOf(availableShapes_all[9]));
        Assert.AreEqual(-1, shapes.IndexOf(availableShapes_all[10]));
        Assert.AreEqual(-1, shapes.IndexOf(availableShapes_all[11]));
        
        Assert.AreEqual(0, shapes2.IndexOf(availableShapes_all[0]));
        Assert.AreEqual(1, shapes2.IndexOf(availableShapes_all[1]));
        Assert.AreEqual(2, shapes2.IndexOf(availableShapes_all[2]));
        Assert.AreEqual(3, shapes2.IndexOf(availableShapes_all[3]));
        Assert.AreEqual(4, shapes2.IndexOf(availableShapes_all[4]));
        Assert.AreEqual(5, shapes2.IndexOf(availableShapes_all[5]));
        Assert.AreEqual(6, shapes2.IndexOf(availableShapes_all[6]));
        Assert.AreEqual(-1, shapes2.IndexOf(availableShapes_all[7]));
        Assert.AreEqual(-1, shapes2.IndexOf(availableShapes_all[8]));
        Assert.AreEqual(7, shapes2.IndexOf(availableShapes_all[9]));
        Assert.AreEqual(8, shapes2.IndexOf(availableShapes_all[10]));
        Assert.AreEqual(9, shapes2.IndexOf(availableShapes_all[11]));
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.rect, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(61, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-200, area);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(10, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.rect, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, .5);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle_rounding.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(false, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(121, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-(Math.PI * 10 * 10) - area), 0.15);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.rect, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rectangle_bias.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(105, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-22*30, area);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(22, bounds.Width);
        Assert.AreEqual(30, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.rect, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "notL.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(61, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-200, area);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(10, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(true, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(69, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-250, area);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_inner.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(68, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-252.8 - area), 0.01);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_ref.shapeIndex);
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.AreEqual(68, clean_ref.Count);
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.LessOrEqual(Math.Abs(-252.8 - area_ref), 0.01);
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.AreEqual(15, bounds_ref.Width);
        Assert.AreEqual(20, bounds_ref.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, iCV:5, iCVariation_scalar:1.0);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_innervariation.svg", FillRule.NonZero, 800, 800, 10);
        // Check point count - start and end points are the same.
        Assert.AreEqual(68, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(area_ref - area), 0.01);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_outer.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(60, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-225.18 - area), 0.01);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_ref.shapeIndex);
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.AreEqual(60, clean_ref.Count);
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.LessOrEqual(Math.Abs(-225.18 - area_ref), 0.01);
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.AreEqual(15, bounds_ref.Width);
        Assert.AreEqual(20, bounds_ref.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, oCV:5, oCVariation_scalar:1.0);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_outervariation.svg", FillRule.NonZero, 800, 800, 10);
        // Check point count - start and end points are the same.
        Assert.AreEqual(60, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(area_ref - area), 0.01);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(59, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-228 - area), 0.015);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_ref.shapeIndex);
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.AreEqual(59, clean_ref.Count);
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.LessOrEqual(Math.Abs(-228 - area_ref), 0.015);
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.AreEqual(15, bounds_ref.Width);
        Assert.AreEqual(20, bounds_ref.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1, oCV:5, oCVariation_scalar:1.0, iCV:10, iCVariation_scalar:1.0);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, out_ref);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_inneroutervar.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(59, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(area_ref - area), 0.015);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_0.shapeIndex);
        PathD out_0 = shape_0.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_0 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension0.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_0 = GeoWrangler.removeDuplicates(out_0);
        // Check point count - start and end points are the same.
        Assert.AreEqual(59, clean_0.Count);
        // Check expected area
        double area_0 = Clipper.Area(out_0);
        Assert.LessOrEqual(Math.Abs(-228 - area_0), 0.015);
        RectD bounds = Clipper.GetBounds(out_0);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_05.shapeIndex);
        PathD out_05 = shape_05.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_05 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension05.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_05 = GeoWrangler.removeDuplicates(out_05);
        // Check point count - start and end points are the same.
        Assert.AreEqual(60, clean_05.Count);
        // Check expected area
        double area_05 = Clipper.Area(out_05);
        Assert.LessOrEqual(Math.Abs(-233.78 - area_05), 0.015);
        RectD bounds_05 = Clipper.GetBounds(out_05);
        Assert.AreEqual(15, bounds_05.Width);
        Assert.AreEqual(20, bounds_05.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_10.shapeIndex);
        PathD out_10 = shape_10.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_10 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tension10.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_10 = GeoWrangler.removeDuplicates(out_10);
        // Check point count - start and end points are the same.
        Assert.AreEqual(60, clean_10.Count);
        // Check expected area
        double area_10 = Clipper.Area(out_10);
        Assert.LessOrEqual(Math.Abs(-239.31 - area_10), 0.015);
        RectD bounds_10 = Clipper.GetBounds(out_10);
        Assert.AreEqual(15, bounds_10.Width);
        Assert.AreEqual(20, bounds_10.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-(3600-((20*20)+(10*15))) - area), 0.0001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(259, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-2915.0441 - area), 0.0001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Tshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "t.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(241, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-((20*60) + (40*10)), area);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Tshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "t_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(204, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1503.028 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Ushape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(true, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(259, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-((60*40) - (30*30)), area);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(40, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Ushape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(false, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(225, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1384.031 - area), 0.001);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(40, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Xshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(197, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1400 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(40, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Xshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(164, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1327.185 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(40, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_ref.shapeIndex);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ref }, true);
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_sbias.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(325, clean.Count);
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(636 - Math.Abs(area_ref - area)), 0.0001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(64, bounds.Width);
        Assert.AreEqual(64, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tip.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(true, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(277, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-2360, area);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(66, bounds.Width);
        Assert.AreEqual(72, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_ref.shapeIndex);
        PathD out_ref = shape_ref.processCorners(false, false, 90, 1, 1);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_ref = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_ref), 0.001);
        Assert.AreEqual(true, orthogonal_ref);
        // Corners can have duplicate points.
        PathD clean_ref = GeoWrangler.removeDuplicates(out_ref);
        // Check point count - start and end points are the same.
        Assert.AreEqual(277, clean_ref.Count);
        // Check expected area
        double area_ref = Clipper.Area(out_ref);
        Assert.AreEqual(-2360, area_ref);
        RectD bounds_ref = Clipper.GetBounds(out_ref);
        Assert.AreEqual(66, bounds_ref.Width);
        Assert.AreEqual(72, bounds_ref.Height);
        
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
        shape_htv.computeCage(0,6, 0);
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_htv.shapeIndex);
        PathD out_htv = shape_htv.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_htv }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tipvarh.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_htv = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_htv), 0.001);
        Assert.AreEqual(true, orthogonal_htv);
        // Corners can have duplicate points.
        PathD clean_htv = GeoWrangler.removeDuplicates(out_htv);
        // Check point count - start and end points are the same.
        Assert.AreEqual(277, clean_htv.Count);
        // Check expected area
        double area_htv = Clipper.Area(out_htv);
        Assert.AreEqual(area_ref, area_htv);
        RectD bounds_htv = Clipper.GetBounds(out_htv);
        Assert.AreEqual(66, bounds_htv.Width);
        Assert.AreEqual(72, bounds_htv.Height);
        
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
        shape_vtv.computeCage(12,0, 0);
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Lshape, shape_vtv.shapeIndex);
        PathD out_vtv = shape_vtv.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_vtv }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "l_tipvarv.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal_vtv = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_vtv), 0.001);
        Assert.AreEqual(true, orthogonal_vtv);
        // Corners can have duplicate points.
        PathD clean_vtv = GeoWrangler.removeDuplicates(out_vtv);
        // Check point count - start and end points are the same.
        Assert.AreEqual(277, clean_vtv.Count);
        // Check expected area
        double area_vtv = Clipper.Area(out_vtv);
        Assert.AreEqual(-2360, area_vtv);
        RectD bounds_vtv = Clipper.GetBounds(out_vtv);
        Assert.AreEqual(66, bounds_vtv.Width);
        Assert.AreEqual(72, bounds_vtv.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Ushape, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "u_brtip.svg", FillRule.NonZero, 800, 800, 10);
        // Ortho check
        // double[] angles = GeoWrangler.angles(GeoWrangler.stripCollinear(out_), true);
        bool orthogonal = GeoWrangler.orthogonal(GeoWrangler.stripCollinear(out_), 0.001);
        Assert.AreEqual(true, orthogonal);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(239, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        // Tip configuration reduces width and height
        Assert.AreEqual(-((60*40) - (25*20)), area);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(40, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_1.shapeIndex);
        PathD out_1 = shape_1.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_bltips.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_1 = GeoWrangler.removeDuplicates(out_1);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean_1.Count);
        // Check expected area
        double area_1 = Clipper.Area(out_1);
        // In this this, the BL tip setting moves the lower edge upwards, reducing the notch vertically.
        // The left tip setting has no effect on the horizontal dimension.
        Assert.LessOrEqual(Math.Abs(-(3600-((20*18)+(10*15))) - area_1), 0.0001);
        RectD bounds_1 = Clipper.GetBounds(clean_1);
        Assert.AreEqual(60, bounds_1.Width);
        Assert.AreEqual(60, bounds_1.Height);
        
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_2.shapeIndex);
        PathD out_2 = shape_2.processCorners(false, false, 90, 1, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "s_brtips.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_2 = GeoWrangler.removeDuplicates(out_2);
        // Check point count - start and end points are the same.
        Assert.AreEqual(297, clean_2.Count);
        // Check expected area
        double area_2 = Clipper.Area(out_2);
        // In this this, the BR tip setting moves the lower edge upwards, reducing the notch vertically.
        // The right tip setting reduces the notch width.
        Assert.LessOrEqual(Math.Abs(-(3600-((15*18)+(10*15))) - area_2), 0.0001);
        RectD bounds_2 = Clipper.GetBounds(clean_2);
        Assert.AreEqual(60, bounds_2.Width);
        Assert.AreEqual(60, bounds_2.Height);

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
        
        Assert.AreEqual(10, offset.x);
        // Inverse due to Clipper Y axis inversion
        Assert.AreEqual(-20, offset.y);
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
        Assert.AreEqual(0, offset_0.x);
        Assert.AreEqual(0, offset_0.y);
        
        ShapeLibrary shape_0 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_0.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_0.shapeIndex);
        PathD out_0 = shape_0.processCorners(false, false, 90, 1, 1);
        out_0 = GeoWrangler.move(out_0, offset_0.x, -offset_0.y);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_0 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape0_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_0 = GeoWrangler.removeDuplicates(out_0);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean_0.Count);
        // Check expected area
        double area_0 = Clipper.Area(out_0);
        Assert.LessOrEqual(Math.Abs(-(3600-((20*20)+(10*15))) - area_0), 0.0001);
        RectD bounds_0 = Clipper.GetBounds(clean_0);
        Assert.AreEqual(0, bounds_0.left);
        Assert.AreEqual(0, bounds_0.top);
        
        // Change our reference to subshape 1
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeRefIndex, 1);
        PointD offset_1 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.AreEqual(0, offset_1.x);
        Assert.AreEqual(7, offset_1.y);

        ShapeLibrary shape_1 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_1.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_1.shapeIndex);
        PathD out_1 = shape_1.processCorners(false, false, 90, 1, 1);
        out_1 = GeoWrangler.move(out_1, offset_1.x, -offset_1.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape1_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_1 = GeoWrangler.removeDuplicates(out_1);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean_1.Count);
        // Check expected area
        double area_1 = Clipper.Area(out_1);
        Assert.LessOrEqual(Math.Abs(-(3600-((20*20)+(10*15))) - area_1), 0.0001);
        RectD bounds_1 = Clipper.GetBounds(clean_1);
        Assert.AreEqual(0, bounds_1.left);
        Assert.AreEqual(-7, bounds_1.top);

        // Change our reference to subshape 2
        shapeSettings.setInt(ShapeSettings.properties_i.subShapeRefIndex, 2);
        PointD offset_2 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.AreEqual(-45, offset_2.x);
        Assert.AreEqual(41, offset_2.y);

        ShapeLibrary shape_2 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_2.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_2.shapeIndex);
        PathD out_2 = shape_2.processCorners(false, false, 90, 1, 1);
        out_2 = GeoWrangler.move(out_2, offset_2.x, -offset_2.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape2_pos_bl.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_2 = GeoWrangler.removeDuplicates(out_2);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean_2.Count);
        // Check expected area
        double area_2 = Clipper.Area(out_2);
        Assert.LessOrEqual(Math.Abs(-(3600-((20*20)+(10*15))) - area_2), 0.0001);
        RectD bounds_2 = Clipper.GetBounds(clean_2);
        Assert.AreEqual(-45, bounds_2.left);
        Assert.AreEqual(-41, bounds_2.top);
        
        // Change our reference to subshape 2
        shapeSettings.setInt(ShapeSettings.properties_i.posInSubShapeIndex, (int)ShapeSettings.subShapeLocations.TR);
        PointD offset_3 = shapeOffsets.doOffsets(0, shapeSettings);
        Assert.AreEqual(-60, offset_3.x);
        Assert.AreEqual(51, offset_3.y);

        ShapeLibrary shape_3 = new ShapeLibrary(shapeTable, shapeSettings);
        shape_3.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.Sshape, shape_3.shapeIndex);
        PathD out_3 = shape_3.processCorners(false, false, 90, 1, 1);
        out_3 = GeoWrangler.move(out_3, offset_3.x, -offset_3.y);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { out_3 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "subshape2_pos_tr.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean_3 = GeoWrangler.removeDuplicates(out_3);
        // Check point count - start and end points are the same.
        Assert.AreEqual(309, clean_3.Count);
        // Check expected area
        double area_3 = Clipper.Area(out_3);
        Assert.LessOrEqual(Math.Abs(-(3600-((20*20)+(10*15))) - area_3), 0.0001);
        RectD bounds_3 = Clipper.GetBounds(clean_3);
        Assert.AreEqual(-60, bounds_3.left);
        Assert.AreEqual(-51, bounds_3.top);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.GEOCORE, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(61, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.AreEqual(-200, area);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(10, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.GEOCORE, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, .5);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_outer.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        Assert.AreEqual(121, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-(Math.PI * 10 * 10) - area), 0.15);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(20, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.GEOCORE, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_inner.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(68, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-252.8 - area), 0.01);
        RectD bounds = Clipper.GetBounds(out_);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(20, bounds.Height);
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
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.GEOCORE, shape.shapeIndex);
        PathD out_ = shape.processCorners(false, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customortho_tip.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(83, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-((10+5) * (20+7)) - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(15, bounds.Width);
        Assert.AreEqual(27, bounds.Height);
    }
}