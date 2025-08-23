using Clipper2Lib;
using geoAnalysis;
using geoWrangler;

namespace UnitTests;

public class GeoAnalysisTests
{
    private static string root_loc = "/d/development/DesignLibs_GPL/shapeanalysis_out/";

    // [SetUp]
    public static void GeoAnalysisSetup()
    {
        testAreaHandler();
        testAreaHandler_perPoly();
        testDistanceHandler();
        testDistanceHandler_multi();
        testChordHandler();
        testAngleHandler();
    }

    [Test]
    public static void testAreaHandler()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 100,
                100, 100,
                100, 0
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -10, -10,
                -10, 110,
                110, 110,
                110, -10,
            })
        };
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, aH.listOfOutputPoints, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "area.svg", FillRule.NonZero, 800, 800, 10);

        Assert.That(aH.area, Is.EqualTo(100 * 100));
    }

    [Test]
    public static void testAreaHandler_perPoly()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 10,
                10, 10,
                10, 0
            }),
            Clipper.MakePath(new double[]
            {
                20, 20,
                20, 50,
                50, 50,
                50, 20
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -10, -10,
                -10, 110,
                110, 110,
                110, -10,
            })
        };
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, true);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, aH.listOfOutputPoints, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "area_perpoly.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(aH.area, Is.EqualTo(10 * 10));
    }

    [Test]
    public static void testDistanceHandler()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 10,
                10, 10,
                10, 0
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                20, 20,
                20, 50,
                50, 50,
                50, 20
            })
        };
        DistanceHandler dH = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacing, false);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance1.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH.distanceString, Is.EqualTo("14.142135623730951"));
        
        DistanceHandler dH2 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacingOld, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH2.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH2.distanceString, Is.EqualTo("14.142135623730951"));
        
        DistanceHandler dH3 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosure, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH3.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance3.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH3.distanceString, Is.EqualTo("-14.142135623730951"));

        DistanceHandler dH4 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosureOld, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH4.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance4.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH4.distanceString, Is.EqualTo("-14.142135623730951"));
    }
    
    [Test]
    public static void testDistanceHandler_multi()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 10,
                10, 10,
                10, 0
            }),
            Clipper.MakePath(new double[]
            {
                20, 20,
                20, 107,
                105, 107,
                105, 20
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -10, -10,
                -10, 110,
                110, 110,
                110, -10,
            })
        };

        bool overlap = GeoWrangler.anyPartialOverlap(aPaths, bPaths);
        bool overlap2 = GeoWrangler.anyPartialOverlap(bPaths, aPaths);
        
        DistanceHandler dH = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacing, false);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance1_multi.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH.distanceString, Is.EqualTo("-5.830951894845301"));
        
        DistanceHandler dH2 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacingOld, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH2.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance2_multi.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH2.distanceString, Is.EqualTo("-5.830951894845301"));
        
        DistanceHandler dH3 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosure, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH3.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance3_multi.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH3.distanceString, Is.EqualTo("5.830951894845301"));
        
        DistanceHandler dH4 = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosureOld, false);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, aPaths);
        SvgUtils.AddClip(svgSrc, bPaths);
        SvgUtils.AddSolution(svgSrc, dH4.resultPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "distance4_multi.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(dH4.distanceString, Is.EqualTo("5.830951894845301"));
    }
    
    [Test]
    public static void testChordHandler()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -100,-20,
                -100, 20,
                100,20,
                100,-20
            }),
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -30, -50,
                -30, 50,
                30, 50,
                30, -50
            })
        };
        // This only gets the aChordLengths
        ChordHandler cH = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.a);
        Assert.That(cH.aChordLengths[0], Is.EqualTo(60));
        Assert.That(cH.aChordLengths[1], Is.EqualTo(60));
        Assert.That(cH.bChordLengths[0], Is.EqualTo(0));
        Assert.That(cH.bChordLengths[1], Is.EqualTo(0));
        
        // This only gets the bChordLengths
        ChordHandler cH2 = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.b);
        Assert.That(cH2.aChordLengths[0], Is.EqualTo(0));
        Assert.That(cH2.aChordLengths[1], Is.EqualTo(0));
        Assert.That(cH2.bChordLengths[0], Is.EqualTo(40));
        Assert.That(cH2.bChordLengths[1], Is.EqualTo(40));
        
        // This gets the full set
        ChordHandler cH3 = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.b + 1);
        Assert.That(cH3.aChordLengths[0], Is.EqualTo(60));
        Assert.That(cH3.aChordLengths[1], Is.EqualTo(60));
        Assert.That(cH3.bChordLengths[0], Is.EqualTo(40));
        Assert.That(cH3.bChordLengths[1], Is.EqualTo(40));
    }

    [Test]
    public static void testAngleHandler()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 50,
                200,50,
                200, 0
            }),
        };
        PathD bPath = Clipper.MakePath(new double[]
        {
            90, 25,
            90, 75,
            110, 75,
            110, 25
        });
        bPath = geoWrangler.GeoWrangler.Rotate(new(100,50), bPath, 45.0);
        
        PathsD bPaths = new()
        {
            bPath
        };

        angleHandler aH = new angleHandler(aPaths, bPaths);
        Assert.That(Math.Abs(45.0 - aH.minimumIntersectionAngle), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void testAngleHandler_NoIntersection()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 50,
                50, 50,
                50, 0
            }),
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                60, 0,
                60, 50,
                110, 50,
                110, 0
            })
        };

        angleHandler aH = new angleHandler(aPaths, bPaths);
        // Should be 180 degrees when no intersection
        Assert.That(aH.minimumIntersectionAngle, Is.EqualTo(180.0));
    }

    [Test]
    public static void testAreaHandler_EmptyPaths()
    {
        PathsD aPaths = new();
        PathsD bPaths = new();
        
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);
        Assert.That(aH.area, Is.EqualTo(0));
    }

    [Test]
    public static void testDistanceHandler_OverlappingShapes()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 20,
                20, 20,
                20, 0
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                10, 10,
                10, 30,
                30, 30,
                30, 10
            })
        };
        
        DistanceHandler dH = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacing, false);
        // Overlapping shapes should give negative distance
        Assert.That(double.Parse(dH.distanceString), Is.LessThan(0));
    }

    [Test]
    public static void testChordHandler_EdgeCases()
    {
        // Test with very thin shapes
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 1,
                100, 1,
                100, 0
            }),
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                -10, -5,
                -10, 5,
                110, 5,
                110, -5
            })
        };
        
        // Test all calculation modes
        ChordHandler cH_a = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.a);
        ChordHandler cH_b = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.b);
        ChordHandler cH_both = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.b + 1);
        
        Assert.That(cH_a.aChordLengths[0], Is.GreaterThan(0));
        Assert.That(cH_b.bChordLengths[0], Is.GreaterThan(0));
        Assert.That(cH_both.aChordLengths[0], Is.GreaterThan(0));
        Assert.That(cH_both.bChordLengths[0], Is.GreaterThan(0));
    }

    [Test]
    public static void testDistanceHandler_AllModes()
    {
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 10,
                10, 10,
                10, 0
            })
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                15, 0,
                15, 10,
                25, 10,
                25, 0
            })
        };
        
        // Test all distance calculation modes
        DistanceHandler dH_spacing = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacing, false);
        DistanceHandler dH_spacingOld = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.spacingOld, false);
        DistanceHandler dH_enclosure = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosure, false);
        DistanceHandler dH_enclosureOld = new DistanceHandler(false, aPaths, bPaths, (int)DistanceHandler.spacingCalcModes.enclosureOld, false);
        
        // All should give same result for non-overlapping shapes
        Assert.That(dH_spacing.distanceString, Is.EqualTo(dH_spacingOld.distanceString));
        Assert.That(dH_enclosure.distanceString, Is.EqualTo(dH_enclosureOld.distanceString));
        // Distance should be 5 (gap between shapes)
        Assert.That(Math.Abs(5.0 - double.Parse(dH_spacing.distanceString)), Is.LessThanOrEqualTo(0.001));
    }
}