using System;
using System.IO;
using System.Linq;
using Clipper2Lib;
using geoCoreLib;
using geoWrangler;
using NUnit.Framework;
using shapeEngine;

namespace UnitTests;

public class GeocoreShortEdgeTests
{
    private static PathD Stair()
    {
        GeoCore geo = new();
        geo.updateGeoCore(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "LayoutRegression", "stair.GDS"), GeoCore.fileType.gds);
        Assert.That(geo.isValid(), Is.True);
        geo.updateGeometry(geo.activeStructure, 0);
        var paths = geo.points(true);
        Assert.That(paths, Has.Count.EqualTo(1));
        Assert.That(Clipper.GetBounds(paths).Width, Is.EqualTo(350));
        return paths[0];
    }

    private static PathD Contour(PathD path, bool legacy, double shortLength, double transition)
    {
        ShapeSettings settings = new();
        settings.setInt(ShapeSettings.properties_i.enabled, 1);
        settings.setInt(ShapeSettings.properties_i.legacyRounding, legacy ? 1 : 0);
        ShapeLibrary shape = new(Enumerable.Range(0, 12).ToArray(), settings);
        shape.setShape((int)ShapeLibrary.shapeNames_all.GEOCORE, GeoWrangler.clockwiseAndReorderXY(path));
        return shape.processCorners(true, false, 90, 1, 1,
            shortEdgeLength: shortLength, maxShortEdgeLength: transition);
    }

    private static bool HasDiagonal(PathD path) => Enumerable.Range(1, path.Count - 1).Any(i =>
        Math.Abs(path[i].x - path[i - 1].x) > 1e-6 && Math.Abs(path[i].y - path[i - 1].y) > 1e-6);

    [Test]
    public void ActualStair_NanometerThresholdsBelowStepLengths_DoNotCollapse()
    {
        PathD stair = Stair();
        Assert.That(HasDiagonal(Contour(stair, false, 0.05, 2)), Is.False);
    }

    [Test]
    public void ActualStair_ModernContouring_HonorsShortEdgeOptions()
    {
        PathD stair = Stair();
        PathD result = Contour(stair, false, 100, 200);
        PathD expected = contourGen.makeContour(GeoWrangler.clockwiseAndReorderXY(stair), 0, 0, 1, 1, 100, 200, 1, enableParallel: true, edgeTension: 0);
        Assert.That(HasDiagonal(expected), Is.True);
        Assert.That(HasDiagonal(result), Is.True, "GEOCORE must not silently switch to legacy contouring");
        Assert.That(Math.Abs(Clipper.Area(Clipper.Xor(new PathsD { result }, new PathsD { expected }, FillRule.EvenOdd, 8))),
            Is.LessThan(0.001), "The GEOCORE wrapper must use the modern short-edge contour");
    }

    [Test]
    public void StairHole_ModernContouring_ForwardsShortEdgeOptions()
    {
        PathD stair = GeoWrangler.clockwiseAndReorderXY(Stair());
        PathD outer = Clipper.MakePath(new double[] { -100,-100, 500,-100, 500,500, -100,500 });
        PathsD source = GeoWrangler.makeKeyHole(new PathsD { outer }, new PathsD { stair }, false, true);
        Assert.That(source, Has.Count.EqualTo(1));
        PathD result = Contour(source[0], false, 100, 200);
        PathD unchanged = Contour(source[0], false, 0, 0);
        Assert.That(Math.Abs(Clipper.Area(Clipper.Xor(new PathsD { result }, new PathsD { unchanged }, FillRule.EvenOdd, 4))),
            Is.GreaterThan(100), "Short-edge settings must reach the decomposed hole contour");
        Assert.That(GeoWrangler.sliverGapRemoval(new PathsD { result }), Has.Count.EqualTo(2));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void BooleanResult_HonorsShortEdgeSettings(bool hole)
    {
        PathsD stair = [Stair()];
        PathsD box = [Clipper.MakePath(new double[] { -100,-100, -100,500, 500,500, 500,-100, -100,-100 })];
        // BOOLEAN previews pass the library Boolean result into GEOCORE
        // contouring. Exercise that boundary without depending on Variance UI.
        PathsD source = GeoWrangler.customBoolean(
            (int)GeoWrangler.LayerFlag.none, hole ? box : stair,
            (int)(hole ? GeoWrangler.LayerFlag.NOT : GeoWrangler.LayerFlag.none), hole ? stair : box,
            (int)GeoWrangler.booleanOperation.AND, 1, 1.03);
        Assert.That(source, Is.Not.Empty);
        PathsD Process(double shortLength, double transition)
        {
            PathsD contoured = new(source.Select(path => Contour(path, false, shortLength, transition)));
            return ProcessOverlaps.processOverlaps(contoured, contoured.Select(_ => false).ToList(),
                1.03, 1).geometry;
        }
        PathsD baseline = Process(0, 0);
        PathsD result = Process(100, 200);
        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(Math.Abs(Clipper.Area(Clipper.Xor(baseline, result, FillRule.EvenOdd, 4))), Is.GreaterThan(100),
            "Short-edge settings must change the BOOLEAN result's contour");
        Assert.That(result.Any(HasDiagonal), Is.True);
        Assert.That(GeoWrangler.sliverGapRemoval(result), Has.Count.EqualTo(hole ? 2 : 1));
    }

    [Test]
    public void ActualStair_ExplicitLegacyMode_RetainsOrthogonalSteps()
    {
        Assert.That(HasDiagonal(Contour(Stair(), true, 100, 200)), Is.False);
    }
}
