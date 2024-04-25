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
        Assert.LessOrEqual(Math.Abs(45.0 - aH.minimumIntersectionAngle), 0.001);
    }
}