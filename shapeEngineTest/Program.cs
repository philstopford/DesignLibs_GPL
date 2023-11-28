using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;
using shapeEngine;

namespace shapeEngineTest;

internal class Program
{
    private static string root_loc = "/d/development/DesignLibs_GPL/shapeengine_out/";

    static int[] shapeTable = new[] {
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
    private static void Main(string[] args)
    {
        rectangleTest();
        rectangleRoundingTest();
        notLTest();
        lTest();
        lInnerRoundingTest();
        lOuterRoundingTest();
        lInnerOuterRoundingTest();
        sTest();
        sRoundingTest();
        tTest();
        tRoundingTest();
        uTest();
        uRoundingTest();
        xTest();
        xRoundingTest();
    }

    private static void rectangleTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.rect);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.horLength, 10.0m, 0);
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.verLength, 20.0m, 0);
        ShapeLibrary shape = new ShapeLibrary(shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex));
        // Check the shape settings are in the shape.
        Assert.AreEqual((int)ShapeLibrary.shapeNames_all.rect, shape.shapeIndex);
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void rectangleRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, .5);
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
    
    // Set for L, but it's a rectangle. We should get the rectangle.
    private static void notLTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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

    private static void lTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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

    private static void lInnerRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void lOuterRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void lInnerOuterRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void sTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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

    private static void sRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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

    private static void tTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void tRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "t_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(206, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1497.005 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(60, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
    }

    private static void uTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void uRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
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
    
    private static void xTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(218, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1400 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(40, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
    }
    
    private static void xRoundingTest()
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
        PathD out_ = shape.processCorners(true, false, true, 0, 0, 0, 0, 0, false, 0, 0, 0, false, 90, 1, 1);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { out_ }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "x_innerouter.svg", FillRule.NonZero, 800, 800, 10);
        // Corners can have duplicate points.
        PathD clean = GeoWrangler.removeDuplicates(out_);
        // Check point count - start and end points are the same.
        Assert.AreEqual(218, clean.Count);
        // Check expected area
        double area = Clipper.Area(out_);
        Assert.LessOrEqual(Math.Abs(-1325.983 - area), 0.001);
        RectD bounds = Clipper.GetBounds(clean);
        Assert.AreEqual(40, bounds.Width);
        Assert.AreEqual(60, bounds.Height);
    }
}