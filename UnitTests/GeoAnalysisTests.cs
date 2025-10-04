using Clipper2Lib;
using geoAnalysis;
using geoWrangler;
using TestHelpers;

namespace UnitTests;

/// <summary>
/// Comprehensive tests for the geoAnalysis library, which provides geometric analysis
/// capabilities including area calculations, distance measurements, chord analysis,
/// and angle computations between polygonal geometries.
/// </summary>
public class GeoAnalysisTests
{
    private static string root_loc = TestOutput.GetPath("shapeanalysis_out");

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

    #region AreaHandler Tests

    /// <summary>
    /// Tests basic AreaHandler functionality with simple overlapping rectangles.
    /// Verifies total area calculation in 'all' mode.
    /// </summary>
    [Test]
    public static void AreaHandler_SimpleOverlap_CalculatesCorrectArea()
    {
        // Arrange: Create two overlapping rectangles
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
                50, 50,
                50, 150,
                150, 150,
                150, 50
            })
        };

        // Act: Calculate intersection area
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);

        // Assert: Intersection should be 50x50 = 2500
        Assert.That(aH.area, Is.EqualTo(2500).Within(0.001));
        Assert.That(aH.listOfOutputPoints, Is.Not.Null);
        Assert.That(aH.listOfOutputPoints.Count, Is.EqualTo(1));
    }

    /// <summary>
    /// Tests AreaHandler with complete overlap (one rectangle inside another).
    /// </summary>
    [Test]
    public static void AreaHandler_CompleteOverlap_CalculatesCorrectArea()
    {
        // Arrange: Smaller rectangle completely inside larger one
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

        // Act
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);

        // Assert: Full area of smaller rectangle
        Assert.That(aH.area, Is.EqualTo(10000).Within(0.001));
    }

    /// <summary>
    /// Tests AreaHandler in per-polygon mode with multiple polygons.
    /// Should return the area of the smallest intersecting polygon.
    /// </summary>
    [Test]
    public static void AreaHandler_PerPolyMode_ReturnsSmallestPolygon()
    {
        // Arrange: Multiple small polygons with different intersection areas
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
                -5, -5,
                -5, 55,
                55, 55,
                55, -5,
            })
        };

        // Act: Use per-polygon mode
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, true);

        // Assert: Should return area of smaller polygon (10x10 = 100)
        Assert.That(aH.area, Is.EqualTo(100).Within(0.001));
    }

    /// <summary>
    /// Tests AreaHandler with no overlap between polygons.
    /// </summary>
    [Test]
    public static void AreaHandler_NoOverlap_ReturnsZeroArea()
    {
        // Arrange: Two separate rectangles with no overlap
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
                20, 30,
                30, 30,
                30, 20
            })
        };

        // Act
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);

        // Assert: No intersection should yield zero area
        Assert.That(aH.area, Is.EqualTo(0).Within(0.001));
        Assert.That(aH.listOfOutputPoints.Count, Is.EqualTo(0));
    }

    /// <summary>
    /// Tests AreaHandler with empty input paths.
    /// </summary>
    [Test]
    public static void AreaHandler_EmptyPaths_ReturnsZeroArea()
    {
        // Arrange: Empty path collections
        PathsD aPaths = new();
        PathsD bPaths = new();

        // Act
        AreaHandler aH = new AreaHandler(aPaths, bPaths, true, false);

        // Assert: Empty inputs should yield zero area
        Assert.That(aH.area, Is.EqualTo(0).Within(0.001));
        Assert.That(aH.listOfOutputPoints.Count, Is.EqualTo(0));
    }

    #endregion

    #region Constants and Enums Tests

    /// <summary>
    /// Tests that ConstantsGA provides expected tolerance value.
    /// </summary>
    [Test]
    public static void ConstantsGA_Tolerance_HasCorrectValue()
    {
        // Assert: Tolerance should be set to expected precision
        Assert.That(ConstantsGA.tolerance, Is.EqualTo(0.001));
    }

    /// <summary>
    /// Tests that Supported.calcModes enum contains expected values.
    /// </summary>
    [Test]
    public static void Supported_CalcModes_ContainsExpectedValues()
    {
        // Assert: All expected calculation modes are available
        Assert.That(Enum.IsDefined(typeof(Supported.calcModes), Supported.calcModes.area), Is.True);
        Assert.That(Enum.IsDefined(typeof(Supported.calcModes), Supported.calcModes.enclosure_spacing_overlap), Is.True);
        Assert.That(Enum.IsDefined(typeof(Supported.calcModes), Supported.calcModes.chord), Is.True);
        Assert.That(Enum.IsDefined(typeof(Supported.calcModes), Supported.calcModes.angle), Is.True);
    }

    /// <summary>
    /// Tests AreaHandler.areaCalcModes enum values.
    /// </summary>
    [Test]
    public static void AreaHandler_areaCalcModes_ContainsExpectedValues()
    {
        // Assert: Area calculation modes are properly defined
        Assert.That(Enum.IsDefined(typeof(AreaHandler.areaCalcModes), AreaHandler.areaCalcModes.all), Is.True);
        Assert.That(Enum.IsDefined(typeof(AreaHandler.areaCalcModes), AreaHandler.areaCalcModes.perpoly), Is.True);
    }

    #endregion

    #region Original Integration Tests (Preserved)

    /// <summary>
    /// Original area handler test - preserved for compatibility.
    /// Tests basic area calculation and SVG output generation.
    /// </summary>

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

    /// <summary>
    /// Original per-polygon area handler test - preserved for compatibility.
    /// </summary>
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

    #endregion

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
        bPath = geoWrangler.GeoWrangler.Rotate(new(100, 50), bPath, 45.0);

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
        // Overlapping shapes should give zero or negative distance  
        Assert.That(double.Parse(dH.distanceString), Is.LessThanOrEqualTo(0));
    }

    [Test]
    public static void testChordHandler_EdgeCases()
    {
        // Test with intersecting rectangles that should produce chords
        PathsD aPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                10, 10,
                10, 90,
                90, 90,
                90, 10
            }),
        };
        PathsD bPaths = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 100,
                100, 100,
                100, 0
            })
        };

        // Test calculation for element a
        ChordHandler cH_a = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.a);
        ChordHandler cH_b = new ChordHandler(aPaths, bPaths, 0.01, (int)ChordHandler.chordCalcElements.b);

        // Verify that chord handlers are created successfully
        Assert.That(cH_a, Is.Not.Null);
        Assert.That(cH_b, Is.Not.Null);
        Assert.That(cH_a.aChordLengths, Is.Not.Null);
        Assert.That(cH_b.bChordLengths, Is.Not.Null);

        // Check that chord length arrays have expected size
        Assert.That(cH_a.aChordLengths.Length, Is.EqualTo(2));
        Assert.That(cH_b.bChordLengths.Length, Is.EqualTo(2));
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