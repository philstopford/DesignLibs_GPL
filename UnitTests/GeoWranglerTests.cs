using Clipper2Lib;
using geoWrangler;

namespace UnitTests;

/// <summary>
/// Comprehensive tests for the geoWrangler library, which provides geometric processing
/// and manipulation utilities including arrays, boolean operations, transformations,
/// and geometric analysis functions.
/// </summary>
public class GeoWranglerTests
{
    private static string root_loc = TestHelpers.TestPath.Get("geowrangler_out") + Path.DirectorySeparatorChar;

    // [SetUp]
    public static void GeoWranglerSetup()
    {
        {
            array_test();
            clean_and_flatten_test();
            close_test();
            customBoolean();
            customBoolean2();
            customBoolean3();
            fragmenter_test();
            from_soup_test();
            inflate_test();
            invert_test();
            meas_angle_test();
            meas_distance_test();
            min_max_test();
            proximity();
            proximity2();
            query_angles_test();
            query_clockwise_test();
            query_enclosed_test();
            query_extents_test();
            query_midpoint_test();
            query_orthogonal_test();
            raycaster_test();
            reorder_test();
            rotate_test();
            skeleton_test();
            strip_collinear_test();
            translate_test();
            unidirectional_bias();
        }
    }

    #region Constants Tests

    /// <summary>
    /// Tests that geoWrangler Constants have expected values.
    /// </summary>
    [Test]
    public static void Constants_Values_AreCorrect()
    {
        // Assert: Verify all constant values are as expected
        Assert.That(Constants.roundingDecimalPrecision, Is.EqualTo(4));
        Assert.That(Constants.tolerance, Is.EqualTo(0.001));
        Assert.That(Constants.scalar_1E2, Is.EqualTo(100));
        Assert.That(Constants.scalar_1E2_inv, Is.EqualTo(0.01));
        Assert.That(Constants.scalar_1E4, Is.EqualTo(10000));
        Assert.That(Constants.scalar_1E4_inv, Is.EqualTo(0.0001));
    }

    /// <summary>
    /// Tests that inverse scaling constants are mathematically correct.
    /// </summary>
    [Test]
    public static void Constants_InverseScaling_IsCorrect()
    {
        // Assert: Inverse constants should be mathematically correct
        Assert.That(Constants.scalar_1E2 * Constants.scalar_1E2_inv, Is.EqualTo(1.0).Within(0.000001));
        Assert.That(Constants.scalar_1E4 * Constants.scalar_1E4_inv, Is.EqualTo(1.0).Within(0.000001));
    }

    #endregion

    #region Array Generation Tests

    /// <summary>
    /// Tests basic array generation with single path and integer counts.
    /// </summary>
    [Test]
    public static void MakeArray_SinglePath_CreatesCorrectArray()
    {
        // Arrange: Create a simple rectangle
        PathD sourcePath = new PathD
        {
            new PointD(0, 0),
            new PointD(10, 0),
            new PointD(10, 10),
            new PointD(0, 10)
        };
        int xCount = 3;
        int yCount = 2;
        double xPitch = 20.0;
        double yPitch = 30.0;

        // Act: Generate array
        PathsD result = GeoWrangler.makeArray(sourcePath, xCount, xPitch, yCount, yPitch);

        // Assert: Should create 3x2 = 6 paths
        Assert.That(result.Count, Is.EqualTo(6));

        // Verify first path (origin)
        Assert.That(result[0][0].x, Is.EqualTo(0).Within(0.001));
        Assert.That(result[0][0].y, Is.EqualTo(0).Within(0.001));

        // Verify translated paths
        Assert.That(result[1][0].x, Is.EqualTo(0).Within(0.001));    // x=0, y=1
        Assert.That(result[1][0].y, Is.EqualTo(30).Within(0.001));

        Assert.That(result[2][0].x, Is.EqualTo(20).Within(0.001));   // x=1, y=0
        Assert.That(result[2][0].y, Is.EqualTo(0).Within(0.001));
    }

    /// <summary>
    /// Tests array generation with multiple source paths.
    /// </summary>
    [Test]
    public static void MakeArray_MultiplePaths_CreatesCorrectArray()
    {
        // Arrange: Create two source paths
        PathsD sourcePaths = new PathsD
        {
            new PathD
            {
                new PointD(0, 0),
                new PointD(5, 0),
                new PointD(5, 5),
                new PointD(0, 5)
            },
            new PathD
            {
                new PointD(10, 10),
                new PointD(15, 10),
                new PointD(15, 15),
                new PointD(10, 15)
            }
        };
        int xCount = 2;
        int yCount = 2;
        double xPitch = 25.0;
        double yPitch = 25.0;

        // Act
        PathsD result = GeoWrangler.makeArray(sourcePaths, xCount, xPitch, yCount, yPitch);

        // Assert: Should create 2x2 array positions × 2 source paths = 8 total paths
        Assert.That(result.Count, Is.EqualTo(8));

        // Each array position should contain both source paths
        Assert.That(result[0][0].x, Is.EqualTo(0).Within(0.001));    // First path, first position
        Assert.That(result[1][0].x, Is.EqualTo(10).Within(0.001));   // Second path, first position
    }

    /// <summary>
    /// Tests array generation with decimal pitch values.
    /// </summary>
    [Test]
    public static void MakeArray_DecimalPitch_HandlesCorrectly()
    {
        // Arrange
        PathD sourcePath = new PathD
        {
            new PointD(0, 0),
            new PointD(1, 1)
        };
        decimal xPitch = 1.5m;
        decimal yPitch = 2.5m;

        // Act
        PathsD result = GeoWrangler.makeArray(sourcePath, 2, xPitch, 2, yPitch);

        // Assert: Should handle decimal pitches correctly
        Assert.That(result.Count, Is.EqualTo(4));
        Assert.That(result[2][0].x, Is.EqualTo(1.5).Within(0.001));  // x=1, y=0 position
        Assert.That(result[1][0].y, Is.EqualTo(2.5).Within(0.001));  // x=0, y=1 position
    }

    /// <summary>
    /// Tests array generation with zero counts.
    /// </summary>
    [Test]
    public static void MakeArray_ZeroCounts_ReturnsEmpty()
    {
        // Arrange
        PathD sourcePath = new PathD { new PointD(0, 0), new PointD(1, 1) };

        // Act & Assert: Zero counts should return empty arrays
        PathsD result1 = GeoWrangler.makeArray(sourcePath, 0, 10.0, 2, 10.0);
        Assert.That(result1.Count, Is.EqualTo(0));

        PathsD result2 = GeoWrangler.makeArray(sourcePath, 2, 10.0, 0, 10.0);
        Assert.That(result2.Count, Is.EqualTo(0));
    }

    #endregion

    #region Boolean Operations Tests

    /// <summary>
    /// Tests GeoWrangler boolean operation enumeration values.
    /// </summary>
    [Test]
    public static void BooleanOperation_Enum_ContainsExpectedValues()
    {
        // Assert: All expected boolean operations are available
        Assert.That(Enum.IsDefined(typeof(GeoWrangler.booleanOperation), GeoWrangler.booleanOperation.AND), Is.True);
        Assert.That(Enum.IsDefined(typeof(GeoWrangler.booleanOperation), GeoWrangler.booleanOperation.OR), Is.True);
    }

    /// <summary>
    /// Tests LayerFlag enumeration values.
    /// </summary>
    [Test]
    public static void LayerFlag_Enum_ContainsExpectedValues()
    {
        // Assert: All expected layer flags are available
        Assert.That(Enum.IsDefined(typeof(GeoWrangler.LayerFlag), GeoWrangler.LayerFlag.none), Is.True);
        Assert.That(Enum.IsDefined(typeof(GeoWrangler.LayerFlag), GeoWrangler.LayerFlag.NOT), Is.True);
    }

    /// <summary>
    /// Tests custom boolean operation with basic input.
    /// </summary>
    [Test]
    public static void CustomBoolean_BasicOperation_ReturnsResult()
    {
        // Arrange: Create two overlapping rectangles
        PathsD firstLayer = new PathsD
        {
            new PathD
            {
                new PointD(0, 0),
                new PointD(20, 0),
                new PointD(20, 20),
                new PointD(0, 20)
            }
        };

        PathsD secondLayer = new PathsD
        {
            new PathD
            {
                new PointD(10, 10),
                new PointD(30, 10),
                new PointD(30, 30),
                new PointD(10, 30)
            }
        };

        // Act: Perform boolean AND operation
        PathsD result = GeoWrangler.customBoolean(
            (int)GeoWrangler.LayerFlag.none, firstLayer,
            (int)GeoWrangler.LayerFlag.none, secondLayer,
            (int)GeoWrangler.booleanOperation.AND,
            0.01, 0.0);

        // Assert: Should return intersection result
        Assert.That(result, Is.Not.Null);
        // Note: Exact result depends on implementation details, 
        // but should have some geometric content for overlapping inputs
    }

    #endregion

    #region Edge Case Tests

    /// <summary>
    /// Tests array generation with single count in one dimension.
    /// </summary>
    [Test]
    public static void MakeArray_SingleCount_WorksCorrectly()
    {
        // Arrange
        PathD sourcePath = new PathD
        {
            new PointD(0, 0),
            new PointD(5, 5)
        };

        // Act: 1xN and Nx1 arrays
        PathsD result1 = GeoWrangler.makeArray(sourcePath, 1, 10.0, 3, 10.0);
        PathsD result2 = GeoWrangler.makeArray(sourcePath, 3, 10.0, 1, 10.0);

        // Assert
        Assert.That(result1.Count, Is.EqualTo(3));  // 1x3 = 3
        Assert.That(result2.Count, Is.EqualTo(3));  // 3x1 = 3
    }

    /// <summary>
    /// Tests array generation with empty source path.
    /// </summary>
    [Test]
    public static void MakeArray_EmptySource_HandlesGracefully()
    {
        // Arrange
        PathD emptyPath = new PathD();

        // Act
        PathsD result = GeoWrangler.makeArray(emptyPath, 2, 10.0, 2, 10.0);

        // Assert: Should create array of empty paths
        Assert.That(result.Count, Is.EqualTo(4));
        Assert.That(result[0].Count, Is.EqualTo(0));
    }

    #endregion

    #region Integration Test Preservation

    /// <summary>
    /// Original array test - preserved for compatibility.
    /// Tests basic array functionality and SVG output generation.
    /// </summary>

    [Test]
    public static void array_test()
    {
        PathD source = new()
        {
            new(0, 0),
            new(0, 30),
            new(10, 30),
            new(10, 0)
        };

        PathsD arrayed = GeoWrangler.makeArray(source, 2, 20m, 3, 10m);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, arrayed, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "arrayed.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(arrayed.Count, Is.EqualTo(6));
        RectD bounds = Clipper.GetBounds(arrayed);
        Assert.That(bounds.left, Is.EqualTo(0));
        Assert.That(bounds.top, Is.EqualTo(0));
        Assert.That(bounds.right, Is.EqualTo(30));
        Assert.That(bounds.bottom, Is.EqualTo(50));

        PathsD arrayed_arrayed = GeoWrangler.makeArray(arrayed, 2, 40m, 3, 60m);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, arrayed_arrayed, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "arrayed_arrayed.svg", FillRule.NonZero, 800, 800, 10);

        Assert.That(arrayed_arrayed.Count, Is.EqualTo(6 * 6));
        RectD bounds_arrayed = Clipper.GetBounds(arrayed_arrayed);
        Assert.That(bounds_arrayed.left, Is.EqualTo(0));
        Assert.That(bounds_arrayed.top, Is.EqualTo(0));
        Assert.That(bounds_arrayed.right, Is.EqualTo(70));
        Assert.That(bounds_arrayed.bottom, Is.EqualTo(170));
    }

    [Test]
    public static void from_soup_test()
    {
        PathsD init = new()
        {
            Clipper.MakePath(new double[]
            {
                50,50,
                50,80,
                80,80,
                80, 50
            }),
            Clipper.MakePath(new double[]
            {
                40,40,
                0,40,
                0,20,
                40, 20
            })
        };

        PathsD fromSoup = GeoWrangler.fromSoup(init);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, fromSoup, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "from_soup.svg", FillRule.NonZero, 800, 800, 10);
        // Assert.AreEqual(3, fromSoup.Count);
        // Assert.AreEqual((80*80) - ((40*20) + (30*30)), Clipper.Area(fromSoup));
    }

    [Test]
    public static void clean_and_flatten_test()
    {
        PathsD init = new()
        {
            Clipper.MakePath(new double[]
            {
                0,0,
                0,80,
                80,80,
                80, 0
            }),
            Clipper.MakePath(new double[]
            {
                40,40,
                0,40,
                0,20,
                40, 20
            })
        };

        PathsD clean_flat = GeoWrangler.clean_and_flatten(init);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, clean_flat, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "clean_flat.svg", FillRule.NonZero, 800, 800, 10);
        //Assert.AreEqual(1, clean_flat.Count);
        //Assert.AreEqual((80*80) - (40*20), Clipper.Area(clean_flat));
    }

    [Test]
    public static void reorder_test()
    {
        PathD source_cw = new()
        {
            new(10, 10),
            new(6, 15),
            new(4, 15),
            new(0, 10),
            new(-15, 6),
            new(-15, 4),
            new(0, 0),
            new(4, -5),
            new(6, -5),
            new PointD(10, 0),
            new PointD(15, 4),
            new PointD(15, 6)
        };
        bool orig_orientation = GeoWrangler.isClockwise(source_cw);

        // Should get minimum X value, then minimum Y, maintaining original orientation.
        PathD xy_reordered = GeoWrangler.reOrderXY(source_cw);
        Assert.That(GeoWrangler.isClockwise(xy_reordered), Is.EqualTo(orig_orientation));
        Assert.That(xy_reordered[0].x, Is.EqualTo(-15));
        Assert.That(xy_reordered[0].y, Is.EqualTo(4));


        // Should get minimum Y value, then minimum X, maintaining original orientation.
        PathD yx_reordered = GeoWrangler.reOrderYX(source_cw);
        Assert.That(GeoWrangler.isClockwise(yx_reordered), Is.EqualTo(orig_orientation));
        Assert.That(yx_reordered[0].x, Is.EqualTo(4));
        Assert.That(yx_reordered[0].y, Is.EqualTo(-5));

        Path64 source_cw_i = Clipper.ScalePath64(source_cw, 1);
        bool orig_orientation_i = GeoWrangler.isClockwise(source_cw_i);

        // Should get minimum X value, then minimum Y, maintaining original orientation.
        Path64 xy_reordered_i = GeoWrangler.reOrderXY(source_cw_i);
        Assert.That(GeoWrangler.isClockwise(xy_reordered_i), Is.EqualTo(orig_orientation_i));
        Assert.That(xy_reordered_i[0].X, Is.EqualTo(-15));
        Assert.That(xy_reordered_i[0].Y, Is.EqualTo(4));


        // Should get minimum Y value, then minimum X, maintaining original orientation.
        Path64 yx_reordered_i = GeoWrangler.reOrderYX(source_cw_i);
        Assert.That(GeoWrangler.isClockwise(yx_reordered_i), Is.EqualTo(orig_orientation));
        Assert.That(yx_reordered_i[0].X, Is.EqualTo(4));
        Assert.That(yx_reordered_i[0].Y, Is.EqualTo(-5));
    }

    [Test]
    public static void close_test()
    {
        PathD source = new()
        {
            new(10, 0),
            new(10, 30),
            new(0, 30),
            new(0, 0),
        };

        PathD closed = GeoWrangler.close(source);
        Assert.That(source.Count, Is.EqualTo(4));
        Assert.That(closed.Count, Is.EqualTo(5));
        Assert.That(closed[^1].x, Is.EqualTo(closed[0].x));
        Assert.That(closed[^1].y, Is.EqualTo(closed[0].y));

        PathD source2 = new()
        {
            new(0, 0),
            new(10, 0),
            new(10, 30),
            new(0, 30),
        };

        PathD closed2 = GeoWrangler.close(source2);
        Assert.That(source2.Count, Is.EqualTo(4));
        Assert.That(closed.Count, Is.EqualTo(5));
        Assert.That(closed2[^1].x, Is.EqualTo(closed2[0].x));
        Assert.That(closed2[^1].y, Is.EqualTo(closed2[0].y));

        source.Add(new(10, 0));
        PathD closed3 = GeoWrangler.close(source);
        Assert.That(source.Count, Is.EqualTo(5));
        Assert.That(closed3.Count, Is.EqualTo(5));
        Assert.That(closed3[^1].x, Is.EqualTo(closed3[0].x));
        Assert.That(closed3[^1].y, Is.EqualTo(closed3[0].y));

        source2.Add(new(0, 0));
        PathD closed4 = GeoWrangler.close(source2);
        Assert.That(source2.Count, Is.EqualTo(5));
        Assert.That(closed4.Count, Is.EqualTo(5));
        Assert.That(closed4[^1].x, Is.EqualTo(closed4[0].x));
        Assert.That(closed4[^1].y, Is.EqualTo(closed4[0].y));
    }

    [Test]
    public static void inflate_test()
    {
        PathD ray = new()
        {
            new(-10, 0),
            new(-10, 10)
        };

        PathD width_1 = GeoWrangler.inflatePath(ray, 100);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, new() { width_1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inflate_width_1.svg", FillRule.NonZero, 800, 800, 10);
        RectD bounds_1 = Clipper.GetBounds(width_1);
        Assert.That(bounds_1.Width, Is.EqualTo(2));
        Assert.That(bounds_1.Height, Is.EqualTo(12));
        Assert.That(bounds_1.top, Is.EqualTo(-1));
        Assert.That(bounds_1.bottom, Is.EqualTo(11));

        PathD resized = GeoWrangler.resize(width_1, 2.0);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { resized }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inflate_resized.svg", FillRule.NonZero, 800, 800, 10);
        RectD bounds_2 = Clipper.GetBounds(resized);
        Assert.That(bounds_2.Width, Is.EqualTo(2 * bounds_1.Width));
        Assert.That(bounds_2.Height, Is.EqualTo(2 * bounds_1.Height));

        PathD ray2 = new()
        {
            new(-10, -10),
            new(10, 10)
        };
        PathD width_2 = GeoWrangler.inflatePath(ray2, 1);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, new() { width_2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inflate_width_2.svg", FillRule.NonZero, 800, 800, 10);
        RectD bounds_3 = Clipper.GetBounds(width_2);
        Assert.That(bounds_3.Width, Is.EqualTo(20.02));
        // Line vertical dimension should not have changed.
        Assert.That(bounds_3.Height, Is.EqualTo(20.02));
    }

    [Test]
    public static void invert_test()
    {
        double unbounded_area = (double)(Int32.MaxValue) * Int32.MaxValue * 4;
        PathD simple_box = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 10,
            10, 10,
            10, 0,
            0, 0
        });

        Fragmenter f = new(1);
        PathD fragmented_box = f.fragmentPath(simple_box);

        PathsD inverted_box = GeoWrangler.invertTone(simple_box, false);
        PathsD inverted_fragmented_box = GeoWrangler.invertTone(fragmented_box, false);
        PathsD inverted_fragmented_box_cl = GeoWrangler.invertTone(fragmented_box, true);
        PathsD inverted_box_tri = GeoWrangler.invertTone(simple_box, false, true);
        PathsD inverted_fragmented_box_tri = GeoWrangler.invertTone(fragmented_box, false, true);
        PathsD inverted_fragmented_box_cl_tri = GeoWrangler.invertTone(fragmented_box, true, true);

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSolution(svgSrc, inverted_box, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_box.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_fragmented_box, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_fragmented_box.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_fragmented_box_cl, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_fragmented_box_cl.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_box_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_box_tri.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_fragmented_box_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_fragmented_box_tri.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_fragmented_box_cl_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_fragmented_box_cl_tri.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(inverted_box.Count, Is.EqualTo(2));
        Assert.That(inverted_box[0][0].x, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_box[0][0].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_box[0][1].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_box[0][1].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_box[0][2].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_box[0][2].y, Is.EqualTo(Int32.MaxValue));

        Assert.That(inverted_box[1][0].x, Is.EqualTo(0));
        Assert.That(inverted_box[1][0].y, Is.EqualTo(0));
        Assert.That(inverted_box[1][1].x, Is.EqualTo(0));
        Assert.That(inverted_box[1][1].y, Is.EqualTo(10));
        Assert.That(inverted_box[1][2].x, Is.EqualTo(10));
        Assert.That(inverted_box[1][2].y, Is.EqualTo(10));

        double ib_area = Clipper.Area(inverted_box);
        Assert.That(ib_area, Is.EqualTo(unbounded_area - (10 * 10)));

        // Stripped collinear points
        Assert.That(inverted_fragmented_box.Count, Is.EqualTo(2));
        Assert.That(inverted_fragmented_box[0][0].x, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box[0][0].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box[0][1].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box[0][1].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box[0][2].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box[0][2].y, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box[1].Count, Is.EqualTo(5));

        double ifb_area = Clipper.Area(inverted_fragmented_box);
        Assert.That(ifb_area, Is.EqualTo(unbounded_area - (10 * 10)));

        // Retained collinear points
        Assert.That(inverted_fragmented_box_cl.Count, Is.EqualTo(2));
        Assert.That(inverted_fragmented_box_cl[0][0].x, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[0][0].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[0][1].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[0][1].y, Is.EqualTo(-Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[0][2].x, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[0][2].y, Is.EqualTo(Int32.MaxValue));
        Assert.That(inverted_fragmented_box_cl[1][0].x, Is.EqualTo(0));
        Assert.That(inverted_fragmented_box_cl[1][0].y, Is.EqualTo(0));
        Assert.That(inverted_fragmented_box_cl[1][1].x, Is.EqualTo(0));
        Assert.That(inverted_fragmented_box_cl[1][1].y, Is.EqualTo(1));
        Assert.That(inverted_fragmented_box_cl[1].Count, Is.EqualTo(41));

        double ifbcl_area = Clipper.Area(inverted_fragmented_box_cl);
        Assert.That(ifbcl_area, Is.EqualTo(unbounded_area - (10 * 10)));

        // Triangulation tests.

        PathsD boxes = new();
        boxes.Add(f.fragmentPath(simple_box));
        boxes.Add(f.fragmentPath(Clipper.MakePath(new double[]
        {
            30, 30,
            30, 40,
            40, 40,
            40, 30,
            30, 30
        })));

        PathsD inverted_boxes = GeoWrangler.invertTone(boxes, false, false);
        PathsD inverted_boxes_cl = GeoWrangler.invertTone(boxes, true, false);
        PathsD inverted_boxes_bounds = GeoWrangler.invertTone(boxes, false, false, true);
        PathsD inverted_boxes_cl_bounds = GeoWrangler.invertTone(boxes, true, false, true);
        PathsD inverted_boxes_tri = GeoWrangler.invertTone(boxes, false, true);
        PathsD inverted_boxes_cl_tri = GeoWrangler.invertTone(boxes, true, true);
        PathsD inverted_boxes_bounds_tri = GeoWrangler.invertTone(boxes, false, true, true);
        PathsD inverted_boxes_cl_bounds_tri = GeoWrangler.invertTone(boxes, true, true, true);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_cl, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_cl.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_bounds, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_bounds.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_cl_bounds, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_cl_bounds.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_tri.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_cl_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_cl_tri.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_bounds_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_bounds_tri.svg", FillRule.EvenOdd, 800, 800, 10);

        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, inverted_boxes_cl_bounds_tri, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "inverted_boxes_cl_bounds_tri.svg", FillRule.EvenOdd, 800, 800, 10);

        // Do we have the area reduction that we expect?
        double area_reduction = 10 * 10 * 2;
        double ibx_area = Clipper.Area(inverted_boxes);
        Assert.That(ibx_area, Is.EqualTo(unbounded_area - area_reduction));
        double ibxcl_area = Clipper.Area(inverted_boxes_cl);
        Assert.That(ibxcl_area, Is.EqualTo(unbounded_area - area_reduction));
        double ibxb_area = Clipper.Area(inverted_boxes_bounds);
        Assert.That(ibxb_area, Is.EqualTo((40 * 40) - area_reduction));
        double ibxclb_area = Clipper.Area(inverted_boxes_cl_bounds);
        Assert.That(ibxclb_area, Is.EqualTo((40 * 40) - area_reduction));
        double ibx_area_tri = Clipper.Area(inverted_boxes_tri);
        Assert.That(ibx_area_tri, Is.EqualTo(-unbounded_area - area_reduction));
        double ibxcl_area_tri = Clipper.Area(inverted_boxes_cl_tri);
        Assert.That(ibxcl_area_tri, Is.EqualTo(-unbounded_area - area_reduction));
        double ibxb_area_tri = Clipper.Area(inverted_boxes_bounds_tri);
        Assert.That(ibxb_area_tri, Is.EqualTo((40 * 40) - area_reduction));
        double ibxclb_area_tri = Clipper.Area(inverted_boxes_cl_bounds_tri);
        Assert.That(ibxclb_area_tri, Is.EqualTo((40 * 40) - area_reduction));
    }

    [Test]
    public static void meas_angle_test()
    {
        PathD _90 = Clipper.MakePath(new double[]
        {
            0, 1,
            0, 0,
            1, 0
        });

        PathD _m90 = Clipper.MakePath(new double[]
        {
            0, 1,
            0, 0,
            -1, 0
        });

        PathD _r90 = Clipper.MakePath(new double[]
        {
            -1, 1,
            0, 0,
            1, 1
        });

        PathD _r45 = Clipper.MakePath(new double[]
        {
            1, 1,
            0, 0,
            1, 0
        });

        PathD _mr45 = Clipper.MakePath(new double[]
        {
            -1, 1,
            0, 0,
            -1, 0
        });

        double angle_1 = GeoWrangler.angleBetweenPoints(_90[0], _90[2], _90[1], true);
        double angle_2 = GeoWrangler.angleBetweenPoints(_m90[0], _m90[2], _m90[1], true);
        double angle_3 = GeoWrangler.angleBetweenPoints(_r90[0], _r90[2], _r90[1], true);
        double angle_4 = GeoWrangler.angleBetweenPoints(_r45[0], _r45[2], _r45[1], true);
        double angle_5 = GeoWrangler.angleBetweenPoints(_mr45[0], _mr45[2], _mr45[1], true);
        Assert.That(Math.Abs(90 - angle_1), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(90 - angle_2), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(90 - angle_3), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45 - angle_4), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45 - angle_5), Is.LessThanOrEqualTo(0.001));

        // Compare with integer handling.
        Path64 _90i = Clipper.ScalePath64(_90, 1.0);
        Path64 _m90i = Clipper.ScalePath64(_m90, 1.0);
        Path64 _r90i = Clipper.ScalePath64(_r90, 1.0);
        Path64 _r45i = Clipper.ScalePath64(_r45, 1.0);
        Path64 _mr45i = Clipper.ScalePath64(_mr45, 1.0);
        double angle_1i = GeoWrangler.angleBetweenPoints(_90i[0], _90i[2], _90i[1], true);
        double angle_2i = GeoWrangler.angleBetweenPoints(_m90i[0], _m90i[2], _m90i[1], true);
        double angle_3i = GeoWrangler.angleBetweenPoints(_r90i[0], _r90i[2], _r90i[1], true);
        double angle_4i = GeoWrangler.angleBetweenPoints(_r45i[0], _r45i[2], _r45i[1], true);
        double angle_5i = GeoWrangler.angleBetweenPoints(_mr45i[0], _mr45i[2], _mr45i[1], true);
        Assert.That(Math.Abs(90 - angle_1i), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(90 - angle_2i), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(90 - angle_3i), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45 - angle_4i), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45 - angle_5i), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void meas_distance_test()
    {
        PointD p1 = new(0, 0);
        PointD p2 = new(0, 1);
        PointD p3 = new(1, 1);
        PointD p4 = new(0, -1);
        PointD p5 = new(1, -1);
        PointD p6 = new(-1, -1);

        double d2 = GeoWrangler.distanceBetweenPoints(p1, p2);
        double d3 = GeoWrangler.distanceBetweenPoints(p1, p3);
        double d4 = GeoWrangler.distanceBetweenPoints(p1, p4);
        double d5 = GeoWrangler.distanceBetweenPoints(p1, p5);
        double d6 = GeoWrangler.distanceBetweenPoints(p1, p6);
        double d7 = GeoWrangler.distanceBetweenPoints(p3, p6);
        PointD d2d = GeoWrangler.PointD_distanceBetweenPoints(p1, p2);
        PointD d3d = GeoWrangler.PointD_distanceBetweenPoints(p1, p3);
        PointD d4d = GeoWrangler.PointD_distanceBetweenPoints(p1, p4);
        PointD d5d = GeoWrangler.PointD_distanceBetweenPoints(p1, p5);
        PointD d6d = GeoWrangler.PointD_distanceBetweenPoints(p1, p6);
        PointD d7d = GeoWrangler.PointD_distanceBetweenPoints(p3, p6);

        double hyp = Math.Sqrt(2);
        Assert.That(d2, Is.EqualTo(1));
        Assert.That(d3, Is.EqualTo(hyp));
        Assert.That(d4, Is.EqualTo(1));
        Assert.That(d5, Is.EqualTo(hyp));
        Assert.That(d6, Is.EqualTo(hyp));
        Assert.That(d7, Is.EqualTo(2 * hyp));

        Assert.That(d2d.x, Is.EqualTo(0));
        Assert.That(d2d.y, Is.EqualTo(-1));
        Assert.That(d3d.x, Is.EqualTo(-1));
        Assert.That(d3d.y, Is.EqualTo(-1));
        Assert.That(d4d.x, Is.EqualTo(0));
        Assert.That(d4d.y, Is.EqualTo(1));
        Assert.That(d5d.x, Is.EqualTo(-1));
        Assert.That(d5d.y, Is.EqualTo(1));
        Assert.That(d6d.x, Is.EqualTo(1));
        Assert.That(d6d.y, Is.EqualTo(1));
        Assert.That(d7d.x, Is.EqualTo(2));
        Assert.That(d7d.y, Is.EqualTo(2));

        // Compare with integer handling.
        Point64 p1i = new(0, 0);
        Point64 p2i = new(0, 1);
        Point64 p3i = new(1, 1);
        Point64 p4i = new(0, -1);
        Point64 p5i = new(1, -1);
        Point64 p6i = new(-1, -1);

        double d2i = GeoWrangler.distanceBetweenPoints(p1i, p2i);
        double d3i = GeoWrangler.distanceBetweenPoints(p1i, p3i);
        double d4i = GeoWrangler.distanceBetweenPoints(p1i, p4i);
        double d5i = GeoWrangler.distanceBetweenPoints(p1i, p5i);
        double d6i = GeoWrangler.distanceBetweenPoints(p1i, p6i);
        double d7i = GeoWrangler.distanceBetweenPoints(p3i, p6i);
        Point64 d2p = GeoWrangler.Point64_distanceBetweenPoints(p1i, p2i);
        Point64 d3p = GeoWrangler.Point64_distanceBetweenPoints(p1i, p3i);
        Point64 d4p = GeoWrangler.Point64_distanceBetweenPoints(p1i, p4i);
        Point64 d5p = GeoWrangler.Point64_distanceBetweenPoints(p1i, p5i);
        Point64 d6p = GeoWrangler.Point64_distanceBetweenPoints(p1i, p6i);
        Point64 d7p = GeoWrangler.Point64_distanceBetweenPoints(p3i, p6i);

        Assert.That(d2i, Is.EqualTo(1));
        Assert.That(d3i, Is.EqualTo(hyp));
        Assert.That(d4i, Is.EqualTo(1));
        Assert.That(d5i, Is.EqualTo(hyp));
        Assert.That(d6i, Is.EqualTo(hyp));
        Assert.That(d7i, Is.EqualTo(2 * hyp));

        Assert.That(d2p.X, Is.EqualTo(-0));
        Assert.That(d2p.Y, Is.EqualTo(-1));
        Assert.That(d3p.X, Is.EqualTo(-1));
        Assert.That(d3p.Y, Is.EqualTo(-1));
        Assert.That(d4p.X, Is.EqualTo(0));
        Assert.That(d4p.Y, Is.EqualTo(1));
        Assert.That(d5p.X, Is.EqualTo(-1));
        Assert.That(d5p.Y, Is.EqualTo(1));
        Assert.That(d6p.X, Is.EqualTo(1));
        Assert.That(d6p.Y, Is.EqualTo(1));
        Assert.That(d7p.X, Is.EqualTo(2));
        Assert.That(d7p.Y, Is.EqualTo(2));
    }

    [Test]
    public static void min_max_test()
    {
        PathD path = Clipper.MakePath(new double[]
        {
            -5, 5,
            5, 10,
            15, 10,
            5, 5
        });

        Path64 pathi = Clipper.ScalePath64(path, 1);

        // Indices of point in path matching query.
        int min_x = GeoWrangler.MinX(path);
        int min_xi = GeoWrangler.MinX(pathi);
        int max_x = GeoWrangler.MaxX(path);
        int max_xi = GeoWrangler.MaxX(pathi);
        int min_y = GeoWrangler.MinY(path);
        int min_yi = GeoWrangler.MinY(pathi);
        int max_y = GeoWrangler.MaxY(path);
        int max_yi = GeoWrangler.MaxY(pathi);
        Assert.That(max_x, Is.EqualTo(2));
        Assert.That(max_xi, Is.EqualTo(2));
        Assert.That(min_x, Is.EqualTo(0));
        Assert.That(min_xi, Is.EqualTo(0));
        Assert.That(max_y, Is.EqualTo(1));
        Assert.That(max_yi, Is.EqualTo(1));
        Assert.That(min_y, Is.EqualTo(0));
        Assert.That(min_yi, Is.EqualTo(0));

        PointD max = GeoWrangler.getMaximumPoint(path);
        Point64 maxi = GeoWrangler.getMaximumPoint(pathi);
        PointD min = GeoWrangler.getMinimumPoint(path);
        Point64 mini = GeoWrangler.getMinimumPoint(pathi);

        Assert.That(max.x, Is.EqualTo(15));
        Assert.That(maxi.X, Is.EqualTo(15));
        Assert.That(min.x, Is.EqualTo(-5));
        Assert.That(mini.X, Is.EqualTo(-5));
        Assert.That(max.y, Is.EqualTo(10));
        Assert.That(maxi.Y, Is.EqualTo(10));
        Assert.That(min.y, Is.EqualTo(5));
        Assert.That(mini.Y, Is.EqualTo(5));
    }

    [Test]
    public static void query_angles_test()
    {
        PathD square = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 10,
            10, 10,
            10, 0,
        });
        PathD square_rev = Clipper.MakePath(new double[]
        {
            0, 0,
            10, 0,
            10, 10,
            0, 10,
        });
        PathD rot_square = GeoWrangler.Rotate(new(5, 5), square, 45);
        PathD rot_square_rev = GeoWrangler.Rotate(new(5, 5), square_rev, 45);

        double[] angles_1 = GeoWrangler.angles(square, true);
        double[] angles_2 = GeoWrangler.angles(square_rev, true);
        double[] angles_3 = GeoWrangler.angles(rot_square, true);
        double[] angles_4 = GeoWrangler.angles(rot_square_rev, true);
        Assert.That(angles_1.Count(), Is.EqualTo(4));
        Assert.That(angles_2.Count(), Is.EqualTo(4));
        Assert.That(angles_3.Count(), Is.EqualTo(4));
        Assert.That(angles_4.Count(), Is.EqualTo(4));
        // Should only have one distinct value, 90, due to pure orthogonal input
        Assert.That(angles_1.Distinct().Count(), Is.EqualTo(1));
        Assert.That(angles_2.Distinct().Count(), Is.EqualTo(1));
        Assert.That(angles_3.Distinct().Count(), Is.EqualTo(1));
        Assert.That(angles_4.Distinct().Count(), Is.EqualTo(1));
        Assert.That(angles_1[0], Is.EqualTo(90));
        Assert.That(angles_2[0], Is.EqualTo(90));
        Assert.That(angles_3[0], Is.EqualTo(90));
        Assert.That(angles_4[0], Is.EqualTo(90));
    }

    [Test]
    public static void query_clockwise_test()
    {
        PathD source = new()
        {
            new(10, 0),
            new(10, 30),
            new(0, 30),
            new(0, 0),
        };

        Assert.That(GeoWrangler.isClockwise(source), Is.False);
        // Clipper's cooordinate system is Y-inverted, so this will be true, annoyingly.
        Assert.That(Clipper.IsPositive(source), Is.True);

        PathD clockwise_path = GeoWrangler.clockwise(source);

        Assert.That(GeoWrangler.isClockwise(clockwise_path), Is.True);
        // Clipper's cooordinate system is Y-inverted, so this will be false, annoyingly.
        Assert.That(Clipper.IsPositive(clockwise_path), Is.False);

        PathD source_cw = new()
        {
            new(10, 10),
            new(6, 15),
            new(4, 15),
            new(0, 10),
            new(-15, 6),
            new(-15, 4),
            new(0, 0),
            new(4, -5),
            new(6, -5),
            new PointD(10, 0),
            new PointD(15, 4),
            new PointD(15, 6)
        };

        // Should get minimum X value, then minimum Y
        PathD xy_cw = GeoWrangler.clockwiseAndReorderXY(source_cw);
        Assert.That(GeoWrangler.isClockwise(xy_cw), Is.True);
        Assert.That(xy_cw[0].x, Is.EqualTo(-15));
        Assert.That(xy_cw[0].y, Is.EqualTo(4));

        // Should get minimum Y value, then minimum X
        PathD yx_cw = GeoWrangler.clockwiseAndReorderYX(source_cw);
        Assert.That(GeoWrangler.isClockwise(yx_cw), Is.True);
        Assert.That(yx_cw[0].x, Is.EqualTo(4));
        Assert.That(yx_cw[0].y, Is.EqualTo(-5));

        PathsD paths = new()
        {
            GeoWrangler.Rotate(new (5, 13), source, 45.0),

            new()
            {
                new(10, 10),
                new(6, 15),
                new(4, 15),
                new(0, 10),
                new(-15, 6),
                new(-15, 4),
                new(0, 0),
                new(4, -5),
                new(6, -5),
                new PointD(10, 0),
                new PointD(15, 4),
                new PointD(15, 6)
            },
        };

        // Should get minimum X value, then minimum Y, clockwise.
        PathsD xy_cw_paths = GeoWrangler.clockwiseAndReorderXY(paths);
        Assert.That(GeoWrangler.isClockwise(xy_cw_paths[0]), Is.True);
        Assert.That(GeoWrangler.isClockwise(xy_cw_paths[1]), Is.True);
        Assert.That(Math.Abs(-10.556 - xy_cw_paths[0][0].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(21.485 - xy_cw_paths[0][0].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(xy_cw_paths[1][0].x, Is.EqualTo(-15));
        Assert.That(xy_cw_paths[1][0].y, Is.EqualTo(4));

        // Should get minimum Y value, then minimum X, clockwise
        PathsD yx_cw_paths = GeoWrangler.clockwiseAndReorderYX(paths);
        Assert.That(GeoWrangler.isClockwise(yx_cw_paths[0]), Is.True);
        Assert.That(GeoWrangler.isClockwise(yx_cw_paths[1]), Is.True);
        Assert.That(Math.Abs(10.656 - yx_cw_paths[0][0].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(0.272 - yx_cw_paths[0][0].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(yx_cw_paths[1][0].x, Is.EqualTo(4));
        Assert.That(yx_cw_paths[1][0].y, Is.EqualTo(-5));
    }

    [Test]
    public static void query_enclosed_test()
    {
        PathD small = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 10,
            10, 10,
            10, 0
        });

        PathD outer = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            20, 20,
            20, 0
        });

        PathD small2 = Clipper.MakePath(new double[]
        {
            15, 15,
            15, 25,
            25, 25,
            25, 15
        });

        PathD small3 = Clipper.MakePath(new double[]
        {
            30, 0,
            30, 10,
            40, 10,
            40, 0
        });

        PathsD large_outer = new()
        {
            Clipper.MakePath(new double[]
            {
                0, 0,
                0, 80,
                80, 80,
                80, 0
            })
        };

        PathsD inners = new()
        {
            small2,
            small3
        };

        PathsD inners2 = new()
        {
            small2,
            Clipper.MakePath(new double[]
            {
                75, 75,
                75, 85,
                85, 85,
                85, 75
            })
        };

        PathsD inners3 = new()
        {
            Clipper.MakePath(new double[]
            {
                -30, 0,
                -30, 10,
                -20, 10,
                -20, 0
            }),
            Clipper.MakePath(new double[]
            {
                85, 10,
                85, 20,
                95, 20,
                95, 10
            })
        };

        PathsD inners4 = new()
        {
            Clipper.MakePath(new double[]
            {
                -30, 0,
                -30, 10,
                -20, 10,
                -20, 0
            }),
            Clipper.MakePath(new double[]
            {
                75, 75,
                75, 85,
                85, 85,
                85, 75
            })
        };

        PathsD inners5 = new()
        {
            Clipper.MakePath(new double[]
            {
                -30, 0,
                -30, 10,
                -20, 10,
                -20, 0
            }),
            Clipper.MakePath(new double[]
            {
                30, 0,
                30, 10,
                40, 10,
                40, 0
            })
        };

        // Full enclosure
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, small);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_small-in-outer.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_small_in_outer = GeoWrangler.anyPartialOverlap(new() { small }, new() { outer });
        // Is small enclosed by outer? Enclosure means no overlap.
        bool enc_small_in_outer = GeoWrangler.enclosed(small, new() { outer }, false);
        bool enc_small_in_outer_strict = GeoWrangler.enclosed(small, new() { outer }, true);

        // Overlap
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, small2);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_small2-in-outer.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_small2_in_outer = GeoWrangler.anyPartialOverlap(new() { small2 }, new() { outer });
        // Is small2 enclosed by outer? Enclosure means no overlap.
        bool enc_small2_in_outer = GeoWrangler.enclosed(small2, new() { outer }, false);
        bool enc_small2_in_outer_strict = GeoWrangler.enclosed(small2, new() { outer }, true);

        // No overlap
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, small3);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_small3-in-outer.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_small3_in_outer = GeoWrangler.anyPartialOverlap(new() { small3 }, new() { outer });
        // Is small3 enclosed by outer? Enclosure means no overlap.
        bool enc_small3_in_outer = GeoWrangler.enclosed(small3, new() { outer }, false);
        bool enc_small3_in_outer_strict = GeoWrangler.enclosed(small3, new() { outer }, true);

        // Reversed queries
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, small);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_outer-in-small.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_outer_in_small = GeoWrangler.anyPartialOverlap(new() { outer }, new() { small });
        // Is outer enclosed by small? Enclosure means no overlap.
        bool enc_outer_in_small = GeoWrangler.enclosed(outer, new() { small }, false);
        bool enc_outer_in_small_strict = GeoWrangler.enclosed(outer, new() { small }, true);

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, small2);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_outer-in-small2.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_outer_in_small2 = GeoWrangler.anyPartialOverlap(new() { outer }, new() { small2 });
        // Is outer enclosed by small2? Enclosure means no overlap.
        bool enc_outer_in_small2 = GeoWrangler.enclosed(outer, new() { small2 }, false);
        bool enc_outer_in_small2_strict = GeoWrangler.enclosed(outer, new() { small2 }, true);

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, small3);
        SvgUtils.AddClip(svgSrc, outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_outer-in-small3.svg", FillRule.NonZero, 800, 800, 10);
        // Overlap? Overlap means no enclosure
        bool olap_outer_in_small3 = GeoWrangler.anyPartialOverlap(new() { outer }, new() { small3 });
        // Is outer enclosed by small3? Enclosure means no overlap.
        bool enc_outer_in_small3 = GeoWrangler.enclosed(outer, new() { small3 }, false);
        bool enc_outer_in_small3_strict = GeoWrangler.enclosed(outer, new() { small3 }, true);

        // Multi-polygon tests
        // Both inners here are fully enclosed in the large outer
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, inners);
        SvgUtils.AddClip(svgSrc, large_outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_inners-in-largeouter.svg", FillRule.NonZero, 800, 800, 10);
        bool enc_largeouter_in_inners = GeoWrangler.enclosed(large_outer, inners, false);
        bool enc_largeouter_in_inners_strict = GeoWrangler.enclosed(large_outer, inners, true);
        bool enc_inners_in_largeouter = GeoWrangler.enclosed(inners, large_outer, false);
        bool enc_inners_in_largeouter_strict = GeoWrangler.enclosed(inners, large_outer, true);
        bool olap_inners_in_largeouter = GeoWrangler.anyPartialOverlap(inners, large_outer);
        bool olap_largeouter_in_inners = GeoWrangler.anyPartialOverlap(large_outer, inners);

        // One inner is fully enclosed, the other has a partial overlap
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, inners2);
        SvgUtils.AddClip(svgSrc, large_outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_inners2-in-largeouter.svg", FillRule.NonZero, 800, 800, 10);
        bool enc_largeouter_in_inners2 = GeoWrangler.enclosed(large_outer, inners2, false);
        bool enc_largeouter_in_inners2_strict = GeoWrangler.enclosed(large_outer, inners2, true);
        bool enc_inners2_in_largeouter = GeoWrangler.enclosed(inners2, large_outer, false);
        bool enc_inners2_in_largeouter_strict = GeoWrangler.enclosed(inners2, large_outer, true);
        bool olap_inners2_in_largeouter = GeoWrangler.anyPartialOverlap(inners2, large_outer);
        bool olap_largeouter_in_inners2 = GeoWrangler.anyPartialOverlap(large_outer, inners2);

        // No inner is within the large outer.
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, inners3);
        SvgUtils.AddClip(svgSrc, large_outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_inners3-in-largeouter.svg", FillRule.NonZero, 800, 800, 10);
        bool enc_largeouter_in_inners3 = GeoWrangler.enclosed(large_outer, inners3, false);
        bool enc_largeouter_in_inners3_strict = GeoWrangler.enclosed(large_outer, inners3, true);
        bool enc_inners3_in_largeouter = GeoWrangler.enclosed(inners3, large_outer, false);
        bool enc_inners3_in_largeouter_strict = GeoWrangler.enclosed(inners3, large_outer, true);
        bool olap_inners3_in_largeouter = GeoWrangler.anyPartialOverlap(inners3, large_outer);
        bool olap_largeouter_in_inners3 = GeoWrangler.anyPartialOverlap(large_outer, inners3);

        // One inner is outside the large outer. One has a partial overlap
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, inners4);
        SvgUtils.AddClip(svgSrc, large_outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_inners4-in-largeouter.svg", FillRule.NonZero, 800, 800, 10);
        bool enc_largeouter_in_inners4 = GeoWrangler.enclosed(large_outer, inners4, false);
        bool enc_largeouter_in_inners4_strict = GeoWrangler.enclosed(large_outer, inners4, true);
        bool enc_inners4_in_largeouter = GeoWrangler.enclosed(inners4, large_outer, false);
        bool enc_inners4_in_largeouter_strict = GeoWrangler.enclosed(inners4, large_outer, true);
        bool olap_inners4_in_largeouter = GeoWrangler.anyPartialOverlap(inners4, large_outer);
        bool olap_largeouter_in_inners4 = GeoWrangler.anyPartialOverlap(large_outer, inners4);

        // One inner is outside the large outer. One has a full enclosure
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, inners5);
        SvgUtils.AddClip(svgSrc, large_outer);
        SvgUtils.SaveToFile(svgSrc, root_loc + "enc_inners5-in-largeouter.svg", FillRule.NonZero, 800, 800, 10);
        bool enc_largeouter_in_inners5 = GeoWrangler.enclosed(large_outer, inners5, false);
        bool enc_largeouter_in_inners5_strict = GeoWrangler.enclosed(large_outer, inners5, true);
        bool enc_inners5_in_largeouter = GeoWrangler.enclosed(inners5, large_outer, false);
        bool enc_inners5_in_largeouter_strict = GeoWrangler.enclosed(inners5, large_outer, true);
        bool olap_inners5_in_largeouter = GeoWrangler.anyPartialOverlap(inners5, large_outer);
        bool olap_largeouter_in_inners5 = GeoWrangler.anyPartialOverlap(large_outer, inners5);

        // Enclosure means no overlap, overlap means no enclosure.
        // As such, these queries should give an opposite signal for any case where there is some form of common area.
        Assert.That(olap_small_in_outer, Is.False);
        Assert.That(enc_small_in_outer, Is.True);
        Assert.That(enc_small_in_outer_strict, Is.True);

        Assert.That(olap_small2_in_outer, Is.True);
        Assert.That(enc_small2_in_outer, Is.False);
        Assert.That(enc_small2_in_outer_strict, Is.False);

        // No common area - thus same response.
        Assert.That(olap_small3_in_outer, Is.False);
        Assert.That(enc_small3_in_outer, Is.False);
        Assert.That(enc_small3_in_outer_strict, Is.False);

        // Reversed queries
        Assert.That(olap_outer_in_small, Is.False);
        Assert.That(enc_outer_in_small, Is.True);
        Assert.That(enc_outer_in_small_strict, Is.True);

        Assert.That(olap_outer_in_small2, Is.True);
        Assert.That(enc_outer_in_small2, Is.False);
        Assert.That(enc_outer_in_small2_strict, Is.False);

        // No common area - thus same response.
        Assert.That(olap_outer_in_small3, Is.False);
        Assert.That(enc_outer_in_small3, Is.False);
        Assert.That(enc_outer_in_small3_strict, Is.False);

        // Multi-polygon case where both inners are fully enclosed.
        Assert.That(olap_inners_in_largeouter, Is.False);
        Assert.That(enc_inners_in_largeouter, Is.True);
        Assert.That(enc_inners_in_largeouter_strict, Is.True);
        // Reverse query
        Assert.That(olap_largeouter_in_inners, Is.False);
        Assert.That(enc_largeouter_in_inners, Is.True);
        Assert.That(enc_largeouter_in_inners_strict, Is.True);

        // Multi-polygon case where there is an overlap for one inner and full enclosure for the other.
        Assert.That(olap_inners2_in_largeouter, Is.True);
        Assert.That(enc_inners2_in_largeouter, Is.False);
        Assert.That(enc_inners2_in_largeouter_strict, Is.False);
        // Reverse query
        Assert.That(olap_largeouter_in_inners2, Is.True);
        Assert.That(enc_largeouter_in_inners2, Is.False);
        Assert.That(enc_largeouter_in_inners2_strict, Is.False);

        // Multi-polygon case where there is no common area
        Assert.That(olap_inners3_in_largeouter, Is.False);
        Assert.That(enc_inners3_in_largeouter, Is.False);
        Assert.That(enc_inners3_in_largeouter_strict, Is.False);
        // Reverse query
        Assert.That(olap_largeouter_in_inners3, Is.False);
        Assert.That(enc_largeouter_in_inners3, Is.False);
        Assert.That(enc_largeouter_in_inners3_strict, Is.False);

        // Multi-polygon case where there is one inner with partial overlap
        Assert.That(olap_inners4_in_largeouter, Is.True);
        Assert.That(enc_inners4_in_largeouter, Is.False);
        Assert.That(enc_inners4_in_largeouter_strict, Is.False);
        // Reverse query
        Assert.That(olap_largeouter_in_inners4, Is.True);
        Assert.That(enc_largeouter_in_inners4, Is.False);
        Assert.That(enc_largeouter_in_inners4_strict, Is.False);

        // Multi-polygon case where there is one inner with full enclosure
        Assert.That(olap_inners5_in_largeouter, Is.True);
        Assert.That(enc_inners5_in_largeouter, Is.False);
        Assert.That(enc_inners5_in_largeouter_strict, Is.False);
        // Reverse query
        Assert.That(olap_largeouter_in_inners5, Is.True);
        Assert.That(enc_largeouter_in_inners5, Is.False);
        Assert.That(enc_largeouter_in_inners5_strict, Is.False);
    }

    [Test]
    public static void query_extents_test()
    {
        PathD test1 = Clipper.MakePath(new double[]
        {
            -15, -30,
            -5, 0,
            -15, 30,
            0, 5,
            15, 30,
            5, 0,
            15, -30,
            0, -5
        });

        PointD test1_extents = GeoWrangler.getExtents(test1);
        Assert.That(test1_extents.x, Is.EqualTo(30));
        Assert.That(test1_extents.y, Is.EqualTo(60));

        PathD test2 = GeoWrangler.move(test1, 30m, 30);

        PointD test2_extents = GeoWrangler.getExtents(test2);
        Assert.That(test2_extents.x, Is.EqualTo(30));
        Assert.That(test2_extents.y, Is.EqualTo(60));

        Path64 test1i = Clipper.ScalePath64(test1, 1);

        Point64 test1i_extents = GeoWrangler.getExtents(test1i);
        Assert.That(test1i_extents.X, Is.EqualTo(30));
        Assert.That(test1i_extents.Y, Is.EqualTo(60));

        Path64 test2i = GeoWrangler.move(test1i, 30, 30);

        Point64 test2i_extents = GeoWrangler.getExtents(test2i);
        Assert.That(test2i_extents.X, Is.EqualTo(30));
        Assert.That(test2i_extents.Y, Is.EqualTo(60));
    }

    [Test]
    public static void query_midpoint_test()
    {
        PathD test1 = Clipper.MakePath(new double[]
        {
            -15, -30,
            -5, 0,
            -15, 30,
            0, 5,
            15, 30,
            5, 0,
            15, -30,
            0, -5
        });

        PointD test1_midpoint = GeoWrangler.midPoint(test1);
        Assert.That(test1_midpoint.x, Is.EqualTo(0));
        Assert.That(test1_midpoint.y, Is.EqualTo(0));

        PathD test2 = GeoWrangler.move(test1, 30m, 60);

        PointD test2_midpoint = GeoWrangler.midPoint(test2);
        Assert.That(test2_midpoint.x, Is.EqualTo(30));
        Assert.That(test2_midpoint.y, Is.EqualTo(60));

        Path64 test1i = Clipper.ScalePath64(test1, 1);

        PointD test1i_midpoint = GeoWrangler.midPoint(test1i);
        Assert.That(test1i_midpoint.x, Is.EqualTo(0));
        Assert.That(test1i_midpoint.y, Is.EqualTo(0));

        Path64 test2i = GeoWrangler.move(test1i, 30, 60);

        PointD test2i_midpoint = GeoWrangler.midPoint(test2i);
        Assert.That(test2i_midpoint.x, Is.EqualTo(30));
        Assert.That(test2i_midpoint.y, Is.EqualTo(60));
    }

    [Test]
    public static void query_orthogonal_test()
    {
        PathD test1 = Clipper.MakePath(new double[]
        {
            -15, -30,
            -5, 0,
            -15, 30,
            0, 5,
            15, 30,
            5, 0,
            15, -30,
            0, -5
        });

        PathD test2 = GeoWrangler.Rotate(GeoWrangler.midPoint(test1), test1, 45);

        Path64 test3 = Clipper.ScalePath64(test1, 1);

        PathD test4 = new(test1);
        test4.Reverse();

        PathD test5 = new(test2);
        test5.Reverse();

        Path64 test6 = new(test3);
        test6.Reverse();

        PathD test7 = Clipper.MakePath(new double[]
        {
            -50, -50,
            -50, 100,
            50, 100,
            50, 0,
            100, 0,
            100, -50,
            0, -50,
            0, -75,
            -25, -75,
            -25, -50
        });

        PathD test8 = GeoWrangler.Rotate(GeoWrangler.midPoint(test7), test7, 45);

        Path64 test9 = Clipper.ScalePath64(test7, 1);

        PathD test10 = new(test7);
        test4.Reverse();

        PathD test11 = new(test8);
        test5.Reverse();

        Path64 test12 = new(test9);
        test6.Reverse();

        Assert.That(GeoWrangler.orthogonal(test1, 0.001), Is.False);
        Assert.That(GeoWrangler.orthogonal(test2, 0.001), Is.False);
        Assert.That(GeoWrangler.orthogonal(test3, 0.001), Is.False);
        Assert.That(GeoWrangler.orthogonal(test4, 0.001), Is.False);
        Assert.That(GeoWrangler.orthogonal(test5, 0.001), Is.False);
        Assert.That(GeoWrangler.orthogonal(test6, 0.001), Is.False);

        Assert.That(GeoWrangler.orthogonal(test7, 0.001), Is.True);
        Assert.That(GeoWrangler.orthogonal(test8, 0.001), Is.True);
        Assert.That(GeoWrangler.orthogonal(test9, 0.001), Is.True);
        Assert.That(GeoWrangler.orthogonal(test10, 0.001), Is.True);
        Assert.That(GeoWrangler.orthogonal(test11, 0.001), Is.True);
        Assert.That(GeoWrangler.orthogonal(test12, 0.001), Is.True);
    }

    [Test]
    public static void rotate_test()
    {
        PathD init = Clipper.MakePath(new double[]
        {
            -15, -10,
            -15,35,
            5, 35,
            5, 0,
            25, 0,
            25, -10
        });

        RectD init_b = Clipper.GetBounds(init);
        Assert.That(init_b.bottom, Is.EqualTo(35));
        Assert.That(init_b.left, Is.EqualTo(-15));
        Assert.That(init_b.Width, Is.EqualTo(40));
        Assert.That(init_b.Height, Is.EqualTo(45));

        PointD init_m = GeoWrangler.midPoint(init);
        Assert.That(init_m.x, Is.EqualTo(5));
        Assert.That(init_m.y, Is.EqualTo(12.5));

        PathD test1 = GeoWrangler.Rotate(init_m, init, 90);
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, new() { test1 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_1.svg", FillRule.NonZero, 800, 800, 10);
        RectD test1_b = Clipper.GetBounds(test1);
        Assert.That(test1_b.bottom, Is.EqualTo(32.5));
        Assert.That(test1_b.left, Is.EqualTo(-17.5));
        Assert.That(test1_b.Width, Is.EqualTo(45));
        Assert.That(test1_b.Height, Is.EqualTo(40));

        PathD test2 = GeoWrangler.Rotate(init[0], init, -90);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, new() { test2 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_2.svg", FillRule.NonZero, 800, 800, 10);
        RectD test2_b = Clipper.GetBounds(test2);
        Assert.That(Math.Abs(-10 - test2_b.bottom), Is.LessThanOrEqualTo(0.001));
        Assert.That(test2_b.left, Is.EqualTo(-15));
        Assert.That(test2_b.Width, Is.EqualTo(45));
        Assert.That(test2_b.Height, Is.EqualTo(40));

        PathD test3 = GeoWrangler.Rotate(init[0], init, -45);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, new() { test3 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_3.svg", FillRule.NonZero, 800, 800, 10);
        RectD test3_b = Clipper.GetBounds(test3);
        Assert.That(Math.Abs(21.819 - test3_b.bottom), Is.LessThanOrEqualTo(0.001));
        Assert.That(test3_b.left, Is.EqualTo(-15));
        Assert.That(Math.Abs(45.961 - test3_b.Width), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(60.104 - test3_b.Height), Is.LessThanOrEqualTo(0.001));

        PathD test4 = GeoWrangler.Rotate(init[2], init, 45);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, new() { test4 }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_4.svg", FillRule.NonZero, 800, 800, 10);
        RectD test4_b = Clipper.GetBounds(test4);
        Assert.That(test4_b.bottom, Is.EqualTo(35));
        Assert.That(Math.Abs(-9.142 - test4_b.left), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(60.104 - test4_b.Width), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45.961 - test4_b.Height), Is.LessThanOrEqualTo(0.001));

        Path64 init_i = Clipper.ScalePath64(init, 1);
        PathD test1_i = GeoWrangler.Rotate(init_m, init_i, 90);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init_i);
        SvgUtils.AddSolution(svgSrc, new() { test1_i }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_1_i.svg", FillRule.NonZero, 800, 800, 10);
        RectD test1_i_b = Clipper.GetBounds(test1_i);
        Assert.That(test1_i_b.bottom, Is.EqualTo(32.5));
        Assert.That(test1_i_b.left, Is.EqualTo(-17.5));
        Assert.That(test1_i_b.Width, Is.EqualTo(45));
        Assert.That(test1_i_b.Height, Is.EqualTo(40));

        PathD test2_i = GeoWrangler.Rotate(new PointD(init_i[0]), init_i, -90);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init_i);
        SvgUtils.AddSolution(svgSrc, new() { test2_i }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_i_2.svg", FillRule.NonZero, 800, 800, 10);
        RectD test2_i_b = Clipper.GetBounds(test2_i);
        Assert.That(Math.Abs(-10 - test2_i_b.bottom), Is.LessThanOrEqualTo(0.001));
        Assert.That(test2_i_b.left, Is.EqualTo(-15));
        Assert.That(test2_i_b.Width, Is.EqualTo(45));
        Assert.That(test2_i_b.Height, Is.EqualTo(40));

        PathD test3_i = GeoWrangler.Rotate(new PointD(init_i[0]), init_i, -45);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init_i);
        SvgUtils.AddSolution(svgSrc, new() { test3_i }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_i_3.svg", FillRule.NonZero, 800, 800, 10);
        RectD test3_i_b = Clipper.GetBounds(test3_i);
        Assert.That(Math.Abs(21.819 - test3_i_b.bottom), Is.LessThanOrEqualTo(0.001));
        Assert.That(test3_i_b.left, Is.EqualTo(-15));
        Assert.That(Math.Abs(45.961 - test3_i_b.Width), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(60.104 - test3_i_b.Height), Is.LessThanOrEqualTo(0.001));

        PathD test4_i = GeoWrangler.Rotate(init[2], init, 45);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, init);
        SvgUtils.AddSolution(svgSrc, new() { test4_i }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rotate_i_4.svg", FillRule.NonZero, 800, 800, 10);
        RectD test4_i_b = Clipper.GetBounds(test4);
        Assert.That(Math.Abs(35 - test4_i_b.bottom), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-9.142 - test4_i_b.left), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(60.104 - test4_i_b.Width), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(45.961 - test4_i_b.Height), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void skeleton_test()
    {
        PathD original = new()
        {
            new(0, 0),
            new(0, 80),
            new(80, 80),
            new(80, 60),
            new(40, 60),
            new(40, 0),
            new(0, 0)
        };

        PathD median = GeoWrangler.skeleton(original);

        SvgWriter svg = new();
        svg.FillRule = FillRule.EvenOdd;

        SvgUtils.AddSubject(svg, original);
        SvgUtils.AddOpenSolution(svg, new PathsD() { median }, true);
        SvgUtils.SaveToFile(svg, root_loc + "median.svg", FillRule.NonZero, 5500, 5500, 10);

        Assert.That(median.Count, Is.EqualTo(6));
        Assert.That(median[0].x, Is.EqualTo(77.5));
        Assert.That(median[0].y, Is.EqualTo(70));
        Assert.That(median[1].x, Is.EqualTo(75));
        Assert.That(median[1].y, Is.EqualTo(70));
        Assert.That(median[2].x, Is.EqualTo(25));
        Assert.That(median[2].y, Is.EqualTo(70));
        Assert.That(median[3].x, Is.EqualTo(20));
        Assert.That(median[3].y, Is.EqualTo(60));
        Assert.That(median[4].x, Is.EqualTo(20));
        Assert.That(median[4].y, Is.EqualTo(13.33));
    }

    [Test]
    public static void unidirectional_bias()
    {
        SvgWriter svg = new();
        svg.FillRule = FillRule.EvenOdd;

        Path64 test = Clipper.MakePath(new[] { 0, 0, 4500, 0, 4500, 4500, 0, 4500 });

        Path64 testvector = Clipper.MakePath(new[] { 0, 0, 500, 0 });

        Paths64 res = Minkowski.Sum(test, testvector, true);

        SvgUtils.AddSubject(svg, test);
        SvgUtils.AddSolution(svg, Clipper.Union(res, new() { test }, FillRule.Positive), true);
        SvgUtils.SaveToFile(svg, root_loc + "unidirectional1.svg", FillRule.NonZero, 5500, 5500, 10);

        Assert.That(Clipper.Area(res), Is.EqualTo(4500000));

        PathD original = new()
        {
            new(0, 0),
            new (0, 50),
            new (20,50),
            new (20, 80),
            new (0,80),
            new(0, 100),
            new(40, 100),
            new(40, 60),
            new(60, 30),
            new(40, 0),
            new(0, 0)
        };

        PathsD original_ = new();
        original_.Add(original);

        PathD vector = new()
        {
            new(0.0, 0.0),
            new PointD(0.0, 10.0)
        };

        PathD vector2 = new()
        {
            new(0.0, 0.0),
            new PointD(0.0, -10.0)
        };

        // Need both transforms to get the full result. Not sure if implementation-related or a bug in C2.
        // No big deal either way.

        PathsD result = Minkowski.Sum(original, vector, true);

        PathsD result_i = Minkowski.Sum(original, vector2, true);

        result.AddRange(result_i);

        SvgUtils.AddSubject(svg, original_);
        SvgUtils.AddSolution(svg, Clipper.Union(result, original_, FillRule.Positive, 2), true);
        SvgUtils.SaveToFile(svg, root_loc + "unidirectional2.svg", FillRule.NonZero, 150, 150, 10);

        Assert.That(Clipper.Area(result), Is.EqualTo(3167));

        PathsD subject = new() { Clipper.Ellipse(new PointD(50.0, 50.0), 20, 20) };

        PathsD result2 = Minkowski.Sum(subject[0], vector, true);

        SvgWriter svg2 = new();
        svg2.FillRule = FillRule.EvenOdd;
        SvgUtils.AddSubject(svg2, subject);
        SvgUtils.AddSolution(svg2, Clipper.Union(result2, subject, FillRule.Positive, 2), true);
        SvgUtils.SaveToFile(svg2, root_loc + "unidirectional3.svg", FillRule.NonZero, 150, 150, 10);

        Assert.That(Math.Abs(Clipper.Area(result2) - 785.5610), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void fragmenter_test()
    {
        PathD original = new()
        {
            new(0, 0),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(0, 0)
        };

        Fragmenter f = new();

        PathD fragmented_10 = f.fragmentPath(original, 1.0);
        PathD fragmented_05 = f.fragmentPath(original, 0.5);

        Assert.That(fragmented_05.Count, Is.EqualTo(81));
        Assert.That(fragmented_10.Count, Is.EqualTo(41));

        // Inject a point into the path that is off-grid for the fragmentation. This has to be collinear to avoid changing the actual shape.
        original.Insert(1, new(0, 5.2));

        // The fragmentation below should change the result on the first edge due to the injected point
        PathD fragmented_b_10 = f.fragmentPath(original, 1.0);
        PathD fragmented_b_05 = f.fragmentPath(original, 0.5);

        Assert.That(fragmented_b_05.Count, Is.EqualTo(80));
        Assert.That(fragmented_b_10.Count, Is.EqualTo(40));

        // The below should give consistent results with the original fragmentation - the injected collinear point should have been stripped.
        PathD refragmented_10 = f.refragmentPath(original, 1.0);
        PathD refragmented_05 = f.refragmentPath(original, 0.5);

        Assert.That(refragmented_05.Count, Is.EqualTo(81));
        Assert.That(refragmented_10.Count, Is.EqualTo(41));
    }

    [Test]
    public static void strip_collinear_test()
    {
        PathD source = new() {
            new(0.02985, 0.18999),
            new(0.00864, 0.21120),
            new(0.01217, 0.21474),
            new(0.01571, 0.21827),
            new(0.02278, 0.21120),
            new(0.02631, 0.21474),
            new(0.02985, 0.21827),
            new(0.00156, 0.24656),
            new(-0.00196, 0.24302),
            new(-0.00550, 0.23949),
            new(-0.01257, 0.24656),
            new(-0.04121, 0.21791),
            new(-0.04829, 0.21084),
            new(-0.05500, 0.20413),
            new(-0.06207, 0.21120),
            new(-0.05500, 0.21827),
            new(-0.06207, 0.22534),
            new(-0.06560, 0.22181),
            new(-0.06914, 0.21828),
            new(-0.07621, 0.22535),
            new(-0.07267, 0.22888),
            new(-0.06914, 0.23242),
            new(-0.06560, 0.23595),
            new(-0.06207, 0.23949),
            new(-0.05853, 0.24302),
            new(-0.05500, 0.24656),
            new(-0.04793, 0.23949),
            new(-0.05499, 0.23242),
            new(-0.05146, 0.22888),
            new(-0.04792, 0.22535),
            new(-0.01964, 0.25363),
            new(-0.02671, 0.26070),
            new(-0.02318, 0.26424),
            new(-0.01964, 0.26777),
            new(-0.06277, 0.31091),
            new(-0.05570, 0.31798),
            new(-0.04792, 0.31020),
            new(-0.04085, 0.31727),
            new(-0.04792, 0.32434),
            new(-0.04085, 0.33141),
            new(-0.01964, 0.31020),
            new(-0.02671, 0.30313),
            new(-0.03378, 0.31020),
            new(-0.04085, 0.30313),
            new(-0.01257, 0.27484),
            new(-0.00550, 0.28191),
            new(0.00156, 0.27484),
            new(0.04399, 0.31727),
            new(0.05106, 0.31020),
            new(0.04753, 0.30666),
            new(0.04399, 0.30313),
            new(0.05106, 0.29606),
            new(0.05813, 0.30313),
            new(0.06520, 0.29606),
            new(0.06167, 0.29252),
            new(0.05813, 0.28899),
            new(0.05460, 0.28545),
            new(0.04399, 0.27484),
            new(0.03692, 0.28191),
            new(0.04399, 0.28899),
            new(0.04045, 0.29252),
            new(0.03692, 0.29606),
            new(0.00864, 0.26777),
            new(0.01571, 0.26070),
            new(0.00864, 0.25363),
            new(0.05177, 0.21050),
            new(0.04470, 0.20343),
            new(0.03692, 0.21120),
            new(0.03338, 0.20767),
            new(0.02985, 0.20413),
            new(0.03692, 0.19706),
            new(0.03692, 0.19706)
        };

        PathD cleaned = GeoWrangler.stripCollinear(source, precision: 6);
        Assert.That(cleaned.Count, Is.EqualTo(71));
    }

    [Test]
    public static void translate_test()
    {
        PathD path = Clipper.MakePath(new double[]
        {
            -5, 5,
            5, 10,
            15, 10,
            5, 5
        });

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, path);

        Path64 pathi = Clipper.ScalePath64(path, 1);

        path = GeoWrangler.move(path, 10m, 20);

        pathi = GeoWrangler.move(pathi, 10, 20);

        SvgUtils.AddSolution(svgSrc, new() { path }, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "translate_path.svg", FillRule.NonZero, 800, 800, 10);

        int min_x = GeoWrangler.MinX(path);
        int min_xi = GeoWrangler.MinX(pathi);
        int max_x = GeoWrangler.MaxX(path);
        int max_xi = GeoWrangler.MaxX(pathi);
        int min_y = GeoWrangler.MinY(path);
        int min_yi = GeoWrangler.MinY(pathi);
        int max_y = GeoWrangler.MaxY(path);
        int max_yi = GeoWrangler.MaxY(pathi);
        Assert.That(max_x, Is.EqualTo(2));
        Assert.That(max_xi, Is.EqualTo(2));
        Assert.That(min_x, Is.EqualTo(0));
        Assert.That(min_xi, Is.EqualTo(0));
        Assert.That(max_y, Is.EqualTo(1));
        Assert.That(max_yi, Is.EqualTo(1));
        Assert.That(min_y, Is.EqualTo(0));
        Assert.That(min_yi, Is.EqualTo(0));

        PointD max = GeoWrangler.getMaximumPoint(path);
        Point64 maxi = GeoWrangler.getMaximumPoint(pathi);
        PointD min = GeoWrangler.getMinimumPoint(path);
        Point64 mini = GeoWrangler.getMinimumPoint(pathi);

        Assert.That(max.x, Is.EqualTo(25));
        Assert.That(maxi.X, Is.EqualTo(25));
        Assert.That(min.x, Is.EqualTo(5));
        Assert.That(mini.X, Is.EqualTo(5));
        Assert.That(max.y, Is.EqualTo(30));
        Assert.That(maxi.Y, Is.EqualTo(30));
        Assert.That(min.y, Is.EqualTo(25));
        Assert.That(mini.Y, Is.EqualTo(25));
    }

    [Test]
    public static void raycaster_test()
    {
        PathD rect_source = Clipper.MakePath(new double[]
        {
            -20, -20,
            -20, 20,
            20, 20,
            20, -20,
            -20, -20
        });

        PathD rect_coll = Clipper.MakePath(new double[]
        {
            -30, -25,
            -30, 45,
            45, 45,
            45, -25,
            -30, -25,
        });

        RayCast rc = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false);
        RayCast rc_project = new RayCast(rect_source, rect_coll, Int32.MaxValue);

        PathsD rc_rays_clipped = rc.getClippedRays();
        PathsD rc_project_rays_clipped = rc_project.getClippedRays();

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_project_rays_clipped.Count, Is.EqualTo(5));
        for (int i = 0; i < rc_project_rays_clipped.Count; i++)
        {
            Assert.That(rc_project_rays_clipped[i][0].x, Is.EqualTo(rect_source[i].x));
            Assert.That(rc_project_rays_clipped[i][0].y, Is.EqualTo(rect_source[i].y));
            Assert.That(rc_project_rays_clipped[i][1].x, Is.EqualTo(0));
            Assert.That(rc_project_rays_clipped[i][1].y, Is.EqualTo(0));
        }

        RayCast rc_nX = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x);
        RayCast rc_project_nX = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.x);

        PathsD rc_nX_rays_clipped = rc_nX.getClippedRays();
        PathsD rc_project_nX_rays_clipped = rc_project_nX.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nX_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nX.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nX_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nX.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nX_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nX_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nX_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_rays_clipped[0][1].y, Is.EqualTo(45));
        Assert.That(rc_nX_rays_clipped[1][1].x, Is.EqualTo(25));
        Assert.That(rc_nX_rays_clipped[1][1].y, Is.EqualTo(-25));
        Assert.That(rc_nX_rays_clipped[2][1].x, Is.EqualTo(-25));
        Assert.That(rc_nX_rays_clipped[2][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nX_rays_clipped[3][1].x), Is.LessThanOrEqualTo(00.001));
        Assert.That(Math.Abs(30 - rc_nX_rays_clipped[3][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nX_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_rays_clipped[4][1].y, Is.EqualTo(45));

        Assert.That(rc_project_nX_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nX_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nX_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nX_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nX_rays_clipped[2][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nX_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nX_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nX_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_nY = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.y);
        RayCast rc_project_nY = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.y);

        PathsD rc_nY_rays_clipped = rc_nY.getClippedRays();
        PathsD rc_project_nY_rays_clipped = rc_project_nY.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nY_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nY.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nY_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nY.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nY_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nY_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nY_rays_clipped[0][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nY_rays_clipped[1][1].x), Is.LessThanOrEqualTo(00.001));
        Assert.That(Math.Abs(30 - rc_nY_rays_clipped[1][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nY_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_nY_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_nY_rays_clipped[3][1].x, Is.EqualTo(25));
        Assert.That(rc_nY_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_nY_rays_clipped[4][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_rays_clipped[4][1].y, Is.EqualTo(-25));

        Assert.That(rc_project_nY_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nY_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nY_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nY_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nY_rays_clipped[2][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nY_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nY_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nY_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_project_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_rays_vertical_clipped = rc_vertical.getClippedRays();
        PathsD rc_project_vertical_rays_clipped = rc_project_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_rays_vertical_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_rays_vertical_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_rays_vertical_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_vertical_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_rays_vertical_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_vertical_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_rays_vertical_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_rays_vertical_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_rays_vertical_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_rays_vertical_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_rays_vertical_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_vertical_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_rays_vertical_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_vertical_clipped[4][1].y, Is.EqualTo(0));

        Assert.That(rc_project_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_vertical_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_vertical_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_vertical_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_vertical_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_vertical_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_vertical_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_vertical_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_vertical_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_vertical_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_vertical_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_vertical_rays_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_project_vertical_rays_clipped[4][1].y, Is.EqualTo(0));

        RayCast rc_nX_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_project_nX_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_nX_vertical_rays_clipped = rc_nX_vertical.getClippedRays();
        PathsD rc_project_nX_vertical_rays_clipped = rc_project_nX_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nX_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nX_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nX_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nX_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nX_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nX_vertical_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_vertical_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_vertical_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_vertical_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_vertical_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_vertical_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_vertical_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_vertical_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_vertical_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_vertical_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nX_vertical_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_vertical_rays_clipped[0][1].y, Is.EqualTo(45));
        Assert.That(rc_nX_vertical_rays_clipped[1][1].x, Is.EqualTo(25));
        Assert.That(rc_nX_vertical_rays_clipped[1][1].y, Is.EqualTo(-25));
        Assert.That(rc_nX_vertical_rays_clipped[2][1].x, Is.EqualTo(-25));
        Assert.That(rc_nX_vertical_rays_clipped[2][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nX_vertical_rays_clipped[3][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(30 - rc_nX_vertical_rays_clipped[3][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nX_vertical_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_vertical_rays_clipped[4][1].y, Is.EqualTo(45));

        Assert.That(rc_project_nX_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nX_vertical_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nX_vertical_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nX_vertical_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_vertical_rays_clipped[1][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nX_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nX_vertical_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_vertical_rays_clipped[3][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nX_vertical_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nX_vertical_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_nY_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_project_nY_vertical = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_nY_vertical_rays_clipped = rc_nY_vertical.getClippedRays();
        PathsD rc_project_nY_vertical_rays_clipped = rc_project_nY_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nY_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nY_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nY_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nY_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nY_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nY_vertical_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_vertical_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_vertical_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_vertical_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_vertical_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_vertical_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_vertical_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_vertical_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_vertical_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_vertical_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nY_vertical_rays_clipped[0][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_vertical_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nY_vertical_rays_clipped[1][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(30 - rc_nY_vertical_rays_clipped[1][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nY_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_nY_vertical_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_nY_vertical_rays_clipped[3][1].x, Is.EqualTo(25));
        Assert.That(rc_nY_vertical_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_nY_vertical_rays_clipped[4][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_vertical_rays_clipped[4][1].y, Is.EqualTo(-25));

        Assert.That(rc_project_nY_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nY_vertical_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nY_vertical_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nY_vertical_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_vertical_rays_clipped[1][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nY_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nY_vertical_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_vertical_rays_clipped[3][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nY_vertical_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nY_vertical_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_project_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_rays_horizontal_clipped = rc_horizontal.getClippedRays();
        PathsD rc_project_horizontal_rays_clipped = rc_project_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_rays_horizontal_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_rays_horizontal_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_rays_horizontal_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_horizontal_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_rays_horizontal_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_horizontal_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_rays_horizontal_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_rays_horizontal_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_rays_horizontal_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_rays_horizontal_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_rays_horizontal_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_rays_horizontal_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_rays_horizontal_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_rays_horizontal_clipped[4][1].y, Is.EqualTo(0));

        Assert.That(rc_project_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_horizontal_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_horizontal_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_horizontal_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_horizontal_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_horizontal_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_horizontal_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_horizontal_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_horizontal_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_horizontal_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_horizontal_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_horizontal_rays_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_project_horizontal_rays_clipped[4][1].y, Is.EqualTo(0));

        RayCast rc_nX_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_project_nX_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_nX_horizontal_rays_clipped = rc_nX_horizontal.getClippedRays();
        PathsD rc_project_nX_horizontal_rays_clipped = rc_project_nX_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nX_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nX_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nX_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nX_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nX_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nX_horizontal_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_horizontal_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_horizontal_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_horizontal_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_horizontal_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_horizontal_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nX_horizontal_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nX_horizontal_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nX_horizontal_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nX_horizontal_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nX_horizontal_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_horizontal_rays_clipped[0][1].y, Is.EqualTo(45));
        Assert.That(rc_nX_horizontal_rays_clipped[1][1].x, Is.EqualTo(25));
        Assert.That(rc_nX_horizontal_rays_clipped[1][1].y, Is.EqualTo(-25));
        Assert.That(rc_nX_horizontal_rays_clipped[2][1].x, Is.EqualTo(-25));
        Assert.That(rc_nX_horizontal_rays_clipped[2][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nX_horizontal_rays_clipped[3][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(30 - rc_nX_horizontal_rays_clipped[3][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nX_horizontal_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_nX_horizontal_rays_clipped[4][1].y, Is.EqualTo(45));

        Assert.That(rc_project_nX_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nX_horizontal_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nX_horizontal_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nX_horizontal_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nX_horizontal_rays_clipped[2][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nX_horizontal_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nX_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nX_horizontal_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nX_horizontal_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_nY_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_project_nY_horizontal = new RayCast(rect_source, rect_coll, Int32.MaxValue, invert: RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_nY_horizontal_rays_clipped = rc_nY_horizontal.getClippedRays();
        PathsD rc_project_nY_horizontal_rays_clipped = rc_project_nY_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_nY_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_nY_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_source);
        SvgUtils.AddClip(svgSrc, rect_coll);
        SvgUtils.AddOpenSolution(svgSrc, rc_project_nY_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_rays_project_nY_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_nY_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_nY_horizontal_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_horizontal_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_horizontal_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_horizontal_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_horizontal_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_horizontal_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_nY_horizontal_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_nY_horizontal_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_nY_horizontal_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_nY_horizontal_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_nY_horizontal_rays_clipped[0][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_horizontal_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(Math.Abs(-30 - rc_nY_horizontal_rays_clipped[1][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(30 - rc_nY_horizontal_rays_clipped[1][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(rc_nY_horizontal_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_nY_horizontal_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_nY_horizontal_rays_clipped[3][1].x, Is.EqualTo(25));
        Assert.That(rc_nY_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_nY_horizontal_rays_clipped[4][1].x, Is.EqualTo(-25));
        Assert.That(rc_nY_horizontal_rays_clipped[4][1].y, Is.EqualTo(-25));

        Assert.That(rc_project_nY_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_project_nY_horizontal_rays_clipped[0][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[0][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[1][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[1][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[2][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[2][0].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[3][0].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[3][0].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[4][0].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[4][0].y, Is.EqualTo(-20));

        Assert.That(rc_project_nY_horizontal_rays_clipped[0][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nY_horizontal_rays_clipped[0][1].y, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[1][1].x, Is.EqualTo(-20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_project_nY_horizontal_rays_clipped[2][1].x, Is.EqualTo(-30));
        Assert.That(rc_project_nY_horizontal_rays_clipped[2][1].y, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[3][1].x, Is.EqualTo(20));
        Assert.That(rc_project_nY_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_project_nY_horizontal_rays_clipped[4][1].x, Is.EqualTo(45));
        Assert.That(rc_project_nY_horizontal_rays_clipped[4][1].y, Is.EqualTo(-20));

        RayCast rc_reversed = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false);
        RayCast rc_reversed_project = new RayCast(rect_coll, rect_source, Int32.MaxValue);

        PathsD rc_reversed_rays_clipped = rc_reversed.getClippedRays();
        PathsD rc_reversed_project_rays_clipped = rc_reversed_project.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_rays_clipped[0][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_rays_clipped[1][1].y, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_rays_clipped[2][1].x, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_rays_clipped[3][1].y, Is.EqualTo(2147483622));
        Assert.That(rc_reversed_project_rays_clipped[4][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nX = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x);
        RayCast rc_reversed_project_nX = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.x);

        PathsD rc_reversed_nX_rays_clipped = rc_reversed_nX.getClippedRays();
        PathsD rc_reversed_project_nX_rays_clipped = rc_reversed_project_nX.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nX_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nX.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nX_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nX.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_nX_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_nX_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_nX_rays_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_rays_clipped[4][1].y, Is.EqualTo(0));

        Assert.That(rc_reversed_project_nX_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_nX_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nX_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nY = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.y);
        RayCast rc_reversed_project_nY = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.y);

        PathsD rc_reversed_nY_rays_clipped = rc_reversed_nY.getClippedRays();
        PathsD rc_reversed_project_nY_rays_clipped = rc_reversed_project_nY.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nY_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nY.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nY_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nY.svg", FillRule.EvenOdd, 800, 800, 10);

        RayCast rc_reversed_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_reversed_project_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_reversed_rays_vertical_clipped = rc_reversed_vertical.getClippedRays();
        PathsD rc_reversed_project_vertical_rays_clipped = rc_reversed_project_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_rays_vertical_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_rays_vertical_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_rays_vertical_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_vertical_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_vertical_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_vertical_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_vertical_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_vertical_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_vertical_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_vertical_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_vertical_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_vertical_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(Math.Abs(-1465267314.6978002 - rc_reversed_rays_vertical_clipped[0][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_vertical_clipped[0][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1465267314.6978002 - rc_reversed_rays_vertical_clipped[1][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1569929278.6048 - rc_reversed_rays_vertical_clipped[1][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1465267329.6978002 - rc_reversed_rays_vertical_clipped[2][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1569929278.6048 - rc_reversed_rays_vertical_clipped[2][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1465267329.6978002 - rc_reversed_rays_vertical_clipped[3][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_vertical_clipped[3][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1465267314.6978002 - rc_reversed_rays_vertical_clipped[4][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_vertical_clipped[4][1].y), Is.LessThanOrEqualTo(0.001));

        Assert.That(rc_reversed_project_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_vertical_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_vertical_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_vertical_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_vertical_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_vertical_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_vertical_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_vertical_rays_clipped[0][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_vertical_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_vertical_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_vertical_rays_clipped[1][1].y, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_vertical_rays_clipped[2][1].x, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_vertical_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_vertical_rays_clipped[3][1].y, Is.EqualTo(2147483622));
        Assert.That(rc_reversed_project_vertical_rays_clipped[4][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_vertical_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nX_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_reversed_project_nX_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_reversed_nX_vertical_rays_clipped = rc_reversed_nX_vertical.getClippedRays();
        PathsD rc_reversed_project_nX_vertical_rays_clipped = rc_reversed_project_nX_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nX_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nX_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nX_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nX_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_nX_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_nX_vertical_rays_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_vertical_rays_clipped[4][1].y, Is.EqualTo(0));

        Assert.That(rc_reversed_project_nX_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_vertical_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nY_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.vertical);
        RayCast rc_reversed_project_nY_vertical = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.vertical);

        PathsD rc_reversed_nY_vertical_rays_clipped = rc_reversed_nY_vertical.getClippedRays();
        PathsD rc_reversed_project_nY_vertical_rays_clipped = rc_reversed_project_nY_vertical.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nY_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nY_vertical.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nY_vertical_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nY_vertical.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_nY_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_nY_vertical_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_vertical_rays_clipped[4][1].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nY_vertical_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_vertical_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_reversed_project_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_reversed_rays_horizontal_clipped = rc_reversed_horizontal.getClippedRays();
        PathsD rc_reversed_project_horizontal_rays_clipped = rc_reversed_project_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_rays_horizontal_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_rays_horizontal_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_rays_horizontal_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_horizontal_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_horizontal_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_horizontal_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_horizontal_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_horizontal_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_horizontal_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_rays_horizontal_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_rays_horizontal_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_rays_horizontal_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(Math.Abs(-1465267314.6978 - rc_reversed_rays_horizontal_clipped[0][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_horizontal_clipped[0][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1465267314.6978 - rc_reversed_rays_horizontal_clipped[1][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1569929278.6048 - rc_reversed_rays_horizontal_clipped[1][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1465267329.6978 - rc_reversed_rays_horizontal_clipped[2][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1569929278.6048 - rc_reversed_rays_horizontal_clipped[2][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(1465267329.6978 - rc_reversed_rays_horizontal_clipped[3][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_horizontal_clipped[3][1].y), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1465267314.6978 - rc_reversed_rays_horizontal_clipped[4][1].x), Is.LessThanOrEqualTo(0.001));
        Assert.That(Math.Abs(-1569929258.6048 - rc_reversed_rays_horizontal_clipped[4][1].y), Is.LessThanOrEqualTo(0.001));

        Assert.That(rc_reversed_project_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_horizontal_rays_clipped[0][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[1][1].y, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[2][1].x, Is.EqualTo(-2147483602));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[3][1].y, Is.EqualTo(2147483622));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[4][1].x, Is.EqualTo(2147483617));
        Assert.That(rc_reversed_project_horizontal_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nX_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_reversed_project_nX_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.x, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_reversed_nX_horizontal_rays_clipped = rc_reversed_nX_horizontal.getClippedRays();
        PathsD rc_reversed_project_nX_horizontal_rays_clipped = rc_reversed_project_nX_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nX_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nX_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nX_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nX_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_nX_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_nX_horizontal_rays_clipped[0][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[0][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[1][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[1][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[2][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[2][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[3][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[3][1].y, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[4][1].x, Is.EqualTo(0));
        Assert.That(rc_reversed_nX_horizontal_rays_clipped[4][1].y, Is.EqualTo(0));

        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nX_horizontal_rays_clipped[4][1].y, Is.EqualTo(-25));

        RayCast rc_reversed_nY_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, projectCorners: false, invert: RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.horizontal);
        RayCast rc_reversed_project_nY_horizontal = new RayCast(rect_coll, rect_source, Int32.MaxValue, invert: RayCast.inversionMode.y, dirOverride: RayCast.forceSingleDirection.horizontal);

        PathsD rc_reversed_nY_horizontal_rays_clipped = rc_reversed_nY_horizontal.getClippedRays();
        PathsD rc_reversed_project_nY_horizontal_rays_clipped = rc_reversed_project_nY_horizontal.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_nY_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_nY_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);
        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, rect_coll);
        SvgUtils.AddClip(svgSrc, rect_source);
        SvgUtils.AddOpenSolution(svgSrc, rc_reversed_project_nY_horizontal_rays_clipped, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "raycaster_reversed_rays_project_nY_horizontal.svg", FillRule.EvenOdd, 800, 800, 10);

        Assert.That(rc_reversed_nY_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_nY_horizontal_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_nY_horizontal_rays_clipped[4][1].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped.Count, Is.EqualTo(5));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[0][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[0][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[1][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[1][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[2][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[2][0].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[3][0].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[3][0].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[4][0].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[4][0].y, Is.EqualTo(-25));

        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[0][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[0][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[1][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[1][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[2][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[2][1].y, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[3][1].x, Is.EqualTo(45));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[3][1].y, Is.EqualTo(-25));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[4][1].x, Is.EqualTo(-30));
        Assert.That(rc_reversed_project_nY_horizontal_rays_clipped[4][1].y, Is.EqualTo(-25));

        PathsD input = new();
        input.Add(Clipper.MakePath(new double[] {
            0, -25,
            0, -24,
            0, -23,
            0, -22,
            0, -21,
            0, -20,
            0, -19,
            0, -18,
            0, -17,
            0, -16,
            0, -15,
            0, -14,
            0, -13,
            0, -12,
            0, -11,
            0, -10,
            0, -9,
            0, -8,
            0, -7,
            0, -6,
            0, -5,
            0.11, -3.96,
            0.43, -2.97,
            0.9500000000000001, -2.06,
            1.6500000000000001, -1.28,
            2.5, -0.67,
            3.45, -0.24,
            4.48, -0.03,
            5, 0,
            6.04, -0.11,
            7.03, -0.43,
            7.94, -0.9500000000000001,
            8.72, -1.6500000000000001,
            9.33, -2.5,
            9.76, -3.45,
            9.97, -4.48,
            10, -5,
            10, -6,
            10, -7,
            10, -8,
            10, -9,
            10, -10,
            10, -11,
            10, -12,
            10, -13,
            10, -14,
            10, -15,
            10, -16,
            10, -17,
            10, -18,
            10, -19,
            10, -20,
            10, -21,
            10, -22,
            10, -23,
            10, -24,
            10, -25,
            9.89, -26.04,
            9.57, -27.03,
            9.05, -27.94,
            8.35, -28.72,
            7.5, -29.330000000000002,
            6.55, -29.76,
            5.5200000000000005, -29.97,
            5, -30,
            3.96, -29.89,
            2.97, -29.57,
            2.06, -29.05,
            1.28, -28.35,
            0.67, -27.5,
            0.24, -26.55,
            0.03, -25.52,
            0, -25,
        }));
        input.Add(Clipper.MakePath(new double[] {
            20, -25,
            20, -24,
            20, -23,
            20, -22,
            20, -21,
            20, -20,
            20, -19,
            20, -18,
            20, -17,
            20, -16,
            20, -15,
            20, -14,
            20, -13,
            20, -12,
            20, -11,
            20, -10,
            20, -9,
            20, -8,
            20, -7,
            20, -6,
            20, -5,
            20.11, -3.96,
            20.43, -2.97,
            20.95, -2.06,
            21.650000000000002, -1.28,
            22.5, -0.67,
            23.45, -0.24,
            24.48, -0.03,
            25, 0,
            26.04, -0.11,
            27.03, -0.43,
            27.94, -0.9500000000000001,
            28.72, -1.6500000000000001,
            29.330000000000002, -2.5,
            29.76, -3.45,
            29.97, -4.48,
            30, -5,
            30, -6,
            30, -7,
            30, -8,
            30, -9,
            30, -10,
            30, -11,
            30, -12,
            30, -13,
            30, -14,
            30, -15,
            30, -16,
            30, -17,
            30, -18,
            30, -19,
            30, -20,
            30, -21,
            30, -22,
            30, -23,
            30, -24,
            30, -25,
            29.89, -26.04,
            29.57, -27.03,
            29.05, -27.94,
            28.35, -28.72,
            27.5, -29.330000000000002,
            26.55, -29.76,
            25.52, -29.97,
            25, -30,
            23.96, -29.89,
            22.97, -29.57,
            22.06, -29.05,
            21.28, -28.35,
            20.67, -27.5,
            20.240000000000002, -26.55,
            20.03, -25.52,
            20, -25,
        }));

        RectD bounds0 = Clipper.GetBounds(input[0]);
        PathD bounds = Clipper.MakePath(new double[]
        {
            bounds0.left, bounds0.bottom,
            bounds0.left, bounds0.top,
            bounds0.right, bounds0.top,
            bounds0.right, bounds0.bottom
        });
        bounds = Clipper.ScalePath(bounds, 2);

        RayCast rc_kidney = new RayCast(input[0], bounds, 128);

        PathsD rc_kidney_rays = rc_kidney.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_rays.svg", FillRule.EvenOdd, 800, 800, 10);

        RayCast rc_kidney_nX = new RayCast(input[0], bounds, 128, invert: RayCast.inversionMode.x);

        PathsD rc_kidney_nX_rays = rc_kidney_nX.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_nX_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_nX_rays.svg", FillRule.EvenOdd, 800, 800, 10);

        RayCast rc_kidney_nY = new RayCast(input[0], bounds, 128, invert: RayCast.inversionMode.y);

        PathsD rc_kidney_nY_rays = rc_kidney_nY.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_nY_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_nY_rays.svg", FillRule.EvenOdd, 800, 800, 10);

        input[0].Reverse();
        RayCast rc_kidney_reversed = new RayCast(input[0], bounds, 128);

        PathsD rc_kidney_reversed_rays = rc_kidney_reversed.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_reversed_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_reversed_rays.svg", FillRule.EvenOdd, 800, 800, 10);

        RayCast rc_kidney_nX_reversed = new RayCast(input[0], bounds, 128, invert: RayCast.inversionMode.x);

        PathsD rc_kidney_nX_reversed_rays = rc_kidney_nX_reversed.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_nX_reversed_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_nX_reversed_rays.svg", FillRule.EvenOdd, 800, 800, 10);

        RayCast rc_kidney_nY_reversed = new RayCast(input[0], bounds, 128, invert: RayCast.inversionMode.y);

        PathsD rc_kidney_nY_reversed_rays = rc_kidney_nY_reversed.getClippedRays();

        svgSrc.ClearAll();
        SvgUtils.AddSubject(svgSrc, input[0]);
        SvgUtils.AddClip(svgSrc, bounds);
        SvgUtils.AddOpenSolution(svgSrc, rc_kidney_nY_reversed_rays, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "rc_kidney_nY_reversed_rays.svg", FillRule.EvenOdd, 800, 800, 10);
    }

    [Test]
    public static void proximity()
    {
        PathsD input = new();
        input.Add(Clipper.MakePath(new[] {
        -5, -15,
        -4.97, -13.950000000000001,
        -4.89, -12.91,
        -4.75, -11.870000000000001,
        -4.5600000000000005, -10.84,
        -4.32, -9.82,
        -4.0200000000000005, -8.82,
        -3.67, -7.83,
        -3.27, -6.87,
        -2.82, -5.92,
        -2.32, -5,
        -1.77, -4.11,
        -1.18, -3.24,
        -0.54, -2.41,
        0.14, -1.62,
        0.86, -0.86,
        1.62, -0.14,
        2.41, 0.54,
        3.24, 1.18,
        4.11, 1.77,
        5, 2.32,
        5.92, 2.82,
        6.87, 3.27,
        7.83, 3.67,
        8.82, 4.0200000000000005,
        9.82, 4.32,
        10.84, 4.5600000000000005,
        11.870000000000001, 4.75,
        12.91, 4.89,
        13.950000000000001, 4.97,
        15, 5,
        16.05, 4.97,
        17.09, 4.89,
        18.13, 4.75,
        19.16, 4.5600000000000005,
        20.18, 4.32,
        21.18, 4.0200000000000005,
        22.17, 3.67,
        23.13, 3.27,
        24.080000000000002, 2.82,
        25, 2.32,
        25.89, 1.77,
        26.76, 1.18,
        27.59, 0.54,
        28.38, -0.14,
        29.14, -0.86,
        29.86, -1.62,
        30.54, -2.41,
        31.18, -3.24,
        31.77, -4.11,
        32.32, -5,
        32.82, -5.92,
        33.27, -6.87,
        33.67, -7.83,
        34.02, -8.82,
        34.32, -9.82,
        34.56, -10.84,
        34.75, -11.870000000000001,
        34.89, -12.91,
        34.97, -13.950000000000001,
        35, -15,
        34.99, -16.05,
        34.96, -17.09,
        34.910000000000004, -18.13,
        34.84, -19.16,
        34.74, -20.18,
        34.63, -21.18,
        34.45, -22.490000000000002,
        34.24, -23.77,
        34, -25,
        33.72, -26.18,
        33.410000000000004, -27.310000000000002,
        33.07, -28.38,
        32.71, -29.39,
        32.32, -30.32,
        31.8, -31.38,
        31.25, -32.32,
        30.55, -33.27,
        29.82, -34.02,
        28.93, -34.63,
        27.89, -34.97,
        27.5, -35,
        26.46, -34.93,
        25.43, -34.71,
        24.45, -34.35,
        23.53, -33.86,
        22.68, -33.25,
        21.93, -32.52,
        21.28, -31.69,
        20.76, -30.79,
        20.37, -29.82,
        20.11, -28.8,
        20, -27.76,
        20, -27.5,
        20.05, -26.46,
        20.19, -25.43,
        20.43, -24.45,
        20.81, -23.42,
        21.28, -22.48,
        21.92, -21.59,
        22.650000000000002, -20.88,
        23.54, -20.330000000000002,
        24.560000000000002, -20.03,
        25, -20,
        26.04, -19.89,
        27.03, -19.57,
        27.94, -19.05,
        28.72, -18.35,
        29.330000000000002, -17.5,
        29.76, -16.55,
        29.97, -15.52,
        30, -15,
        29.89, -13.96,
        29.57, -12.97,
        29.05, -12.06,
        28.35, -11.28,
        27.5, -10.67,
        26.55, -10.24,
        25.52, -10.03,
        25, -10,
        23.96, -9.89,
        22.97, -9.57,
        22.06, -9.05,
        21.28, -8.35,
        20.67, -7.5,
        20.240000000000002, -6.55,
        20.03, -5.5200000000000005,
        20, -5,
        19.89, -3.96,
        19.57, -2.97,
        19.05, -2.06,
        18.35, -1.28,
        17.5, -0.67,
        16.55, -0.24,
        15.52, -0.03,
        15, 0,
        13.96, -0.11,
        12.97, -0.43,
        12.06, -0.9500000000000001,
        11.28, -1.6500000000000001,
        10.67, -2.5,
        10.24, -3.45,
        10.03, -4.48,
        10, -5,
        9.89, -6.04,
        9.57, -7.03,
        9.05, -7.94,
        8.35, -8.72,
        7.5, -9.33,
        6.55, -9.76,
        5.5200000000000005, -9.97,
        5, -10,
        3.96, -10.11,
        2.97, -10.43,
        2.06, -10.950000000000001,
        1.28, -11.65,
        0.67, -12.5,
        0.24, -13.450000000000001,
        0.03, -14.48,
        0, -15,
        0.11, -16.04,
        0.43, -17.03,
        0.9500000000000001, -17.94,
        1.6500000000000001, -18.72,
        2.5, -19.330000000000002,
        3.45, -19.76,
        4.48, -19.97,
        5, -20,
        6.04, -20.16,
        6.95, -20.6,
        7.8, -21.28,
        8.47, -22.1,
        8.99, -22.990000000000002,
        9.41, -23.98,
        9.73, -25.060000000000002,
        9.91, -26.07,
        9.99, -27.11,
        10, -27.5,
        9.93, -28.54,
        9.71, -29.57,
        9.35, -30.55,
        8.86, -31.470000000000002,
        8.25, -32.32,
        7.5200000000000005, -33.07,
        6.69, -33.72,
        5.79, -34.24,
        4.82, -34.63,
        3.8000000000000003, -34.89,
        2.7600000000000002, -35,
        2.5, -35,
        1.46, -34.81,
        0.56, -34.32,
        -0.31, -33.54,
        -1.02, -32.660000000000004,
        -1.58, -31.77,
        -2.12, -30.76,
        -2.61, -29.63,
        -2.99, -28.64,
        -3.33, -27.59,
        -3.64, -26.47,
        -3.93, -25.3,
        -4.18, -24.080000000000002,
        -4.4, -22.81,
        -4.59, -21.51,
        -4.71, -20.51,
        -4.8100000000000005, -19.5,
        -4.89, -18.47,
        -4.94, -17.44,
        -4.98, -16.4,
        -5, -15.35,
        -5, -15,
        }));

        GeometryResult res = Proximity.proximityBias(input, new() { false }, 5, 60, 64, 0, 1, 1.03m, 1, false, 1000);

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, input);
        SvgUtils.AddSolution(svgSrc, res.geometry, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "proximity.svg", FillRule.NonZero, 800, 800, 10);

        Assert.That(res.geometry.Count, Is.EqualTo(1));
        double area = Clipper.Area(res.geometry[0]);
        Assert.That(res.geometry[0].Count, Is.EqualTo(211));
        Assert.That(Math.Abs(area - 1539.759), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void proximity2()
    {
        PathsD input = new();
        input.Add(Clipper.MakePath(new double[] {
            0, -25,
            0, -24,
            0, -23,
            0, -22,
            0, -21,
            0, -20,
            0, -19,
            0, -18,
            0, -17,
            0, -16,
            0, -15,
            0, -14,
            0, -13,
            0, -12,
            0, -11,
            0, -10,
            0, -9,
            0, -8,
            0, -7,
            0, -6,
            0, -5,
            0.11, -3.96,
            0.43, -2.97,
            0.9500000000000001, -2.06,
            1.6500000000000001, -1.28,
            2.5, -0.67,
            3.45, -0.24,
            4.48, -0.03,
            5, 0,
            6.04, -0.11,
            7.03, -0.43,
            7.94, -0.9500000000000001,
            8.72, -1.6500000000000001,
            9.33, -2.5,
            9.76, -3.45,
            9.97, -4.48,
            10, -5,
            10, -6,
            10, -7,
            10, -8,
            10, -9,
            10, -10,
            10, -11,
            10, -12,
            10, -13,
            10, -14,
            10, -15,
            10, -16,
            10, -17,
            10, -18,
            10, -19,
            10, -20,
            10, -21,
            10, -22,
            10, -23,
            10, -24,
            10, -25,
            9.89, -26.04,
            9.57, -27.03,
            9.05, -27.94,
            8.35, -28.72,
            7.5, -29.330000000000002,
            6.55, -29.76,
            5.5200000000000005, -29.97,
            5, -30,
            3.96, -29.89,
            2.97, -29.57,
            2.06, -29.05,
            1.28, -28.35,
            0.67, -27.5,
            0.24, -26.55,
            0.03, -25.52,
            0, -25,
        }));
        input.Add(Clipper.MakePath(new double[] {
            20, -25,
            20, -24,
            20, -23,
            20, -22,
            20, -21,
            20, -20,
            20, -19,
            20, -18,
            20, -17,
            20, -16,
            20, -15,
            20, -14,
            20, -13,
            20, -12,
            20, -11,
            20, -10,
            20, -9,
            20, -8,
            20, -7,
            20, -6,
            20, -5,
            20.11, -3.96,
            20.43, -2.97,
            20.95, -2.06,
            21.650000000000002, -1.28,
            22.5, -0.67,
            23.45, -0.24,
            24.48, -0.03,
            25, 0,
            26.04, -0.11,
            27.03, -0.43,
            27.94, -0.9500000000000001,
            28.72, -1.6500000000000001,
            29.330000000000002, -2.5,
            29.76, -3.45,
            29.97, -4.48,
            30, -5,
            30, -6,
            30, -7,
            30, -8,
            30, -9,
            30, -10,
            30, -11,
            30, -12,
            30, -13,
            30, -14,
            30, -15,
            30, -16,
            30, -17,
            30, -18,
            30, -19,
            30, -20,
            30, -21,
            30, -22,
            30, -23,
            30, -24,
            30, -25,
            29.89, -26.04,
            29.57, -27.03,
            29.05, -27.94,
            28.35, -28.72,
            27.5, -29.330000000000002,
            26.55, -29.76,
            25.52, -29.97,
            25, -30,
            23.96, -29.89,
            22.97, -29.57,
            22.06, -29.05,
            21.28, -28.35,
            20.67, -27.5,
            20.240000000000002, -26.55,
            20.03, -25.52,
            20, -25,
        }));

        GeometryResult res = Proximity.proximityBias(input, new() { false, false }, 5, 30, 64, 0, 1, 1.03m, 1, false, 1000);

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, input);
        SvgUtils.AddSolution(svgSrc, res.geometry, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "proximity2.svg", FillRule.NonZero, 800, 800, 10);

        Assert.That(res.geometry.Count, Is.EqualTo(2));
        double area = Clipper.Area(res.geometry);
        Assert.That(res.geometry[0].Count, Is.EqualTo(73));
        Assert.That(res.geometry[1].Count, Is.EqualTo(73));
        Assert.That(Math.Abs(area - 1179.626), Is.LessThanOrEqualTo(0.001));
    }

    // This tests the complex boolean engine in GeoWrangler, which should deliver a keyholed representation of the 
    // boolean output, if needed.
    // This is a little more sophisticated than just calling Clipper directly - it pulls in various GeoWrangler 
    // systems to work its magic, so this is something of an integrated test.
    // Note that some area is lost due to the keyholed representation.
    [Test]
    public static void customBoolean()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            90,0,
            90,80,
            100,80,
            100,0
        }));
        double aArea = Clipper.Area(layerAPaths);

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));
        double bArea = Clipper.Area(layerBPaths);

        // Subtract layerBPaths from layerAPaths
        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths,
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT,
            secondLayer: layerBPaths,
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not.svg", FillRule.NonZero, 800, 800, 10);

        double notArea_expected = aArea - bArea;
        double notArea = Clipper.Area(booleanPaths);
        Assert.That(booleanPaths.Count, Is.EqualTo(2));
        Assert.That(Math.Abs(notArea_expected - notArea + (2 * (GeoWrangler.keyhole_sizing * 0.001))) - 2, Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void customBoolean2()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            -10,30,
            -10,50,
            90,50,
            90,30
        }));

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));

        // Subtract layerBPaths from layerAPaths
        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths,
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT,
            secondLayer: layerBPaths,
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not2.svg", FillRule.NonZero, 800, 800, 10);

        double notArea_expected = -5598;
        double notArea = Clipper.Area(booleanPaths);
        Assert.That(booleanPaths.Count, Is.EqualTo(1));
        Assert.That(Math.Abs(notArea_expected - notArea + (2 * GeoWrangler.keyhole_sizing * 0.001)), Is.LessThanOrEqualTo(0.001));
    }

    [Test]
    public static void customBoolean3()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double[]
        {
            90,0,
            90,80,
            100,80,
            100,0
        }));
        double aArea = Clipper.Area(layerAPaths);

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));
        double bArea = Clipper.Area(layerBPaths);

        // Subtract layerBPaths from layerAPaths
        PathsD tmpPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths,
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT,
            secondLayer: layerBPaths,
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        PathsD cPaths = new();
        cPaths.Add(Clipper.MakePath(new double[]
        {
            -10,30,
            -10,35,
            90,35,
            90,30
        }));

        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: tmpPaths,
            secondLayerOperator: (int)GeoWrangler.LayerFlag.none,
            secondLayer: cPaths,
            booleanFlag: (int)GeoWrangler.booleanOperation.OR,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not3.svg", FillRule.NonZero, 800, 800, 10);

        double notArea_expected = -6195.989;
        double notArea = Clipper.Area(booleanPaths);
        Assert.That(booleanPaths.Count, Is.EqualTo(1));
        Assert.That(Math.Abs(notArea_expected - notArea), Is.LessThanOrEqualTo(0.0011));
    }

    #endregion
}