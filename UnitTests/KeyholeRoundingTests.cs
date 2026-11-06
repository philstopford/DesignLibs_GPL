using System;
using System.IO;
using System.Linq;
using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;
using shapeEngine;
using System.Text.Json;

namespace UnitTests;

public class KeyholeRoundingTests
{
    [TestCase(1.03)]
    [TestCase(2.06)]
    public void Rekeyholing_ClosesExistingSeamBeforeAddingAnotherHole(double extension)
    {
        PathD outer = Clipper.MakePath(new double[] { 0, 0, 100, 0, 100, 100, 0, 100 });
        PathD firstHole = Clipper.MakePath(new double[] { 40, 40, 40, 60, 60, 60, 60, 40 });
        PathsD original = GeoWrangler.makeKeyHole(new PathsD { outer }, new PathsD { firstHole }, false, true);
        PathD secondHole = Clipper.MakePath(new double[] { 25, 55, 25, 65, 45, 65, 45, 55 });
        PathsD input = new(original) { secondHole };
        PathsD physical = GeoWrangler.sliverGapRemoval(input, customSizing: 2 * GeoWrangler.keyhole_sizing);
        PathsD expected = GeoWrangler.makeKeyHole(physical, false, true, extension: extension);
        PathsD actual = GeoWrangler.makeKeyHole(input, false, true, extension: extension);
        PathsD difference = Clipper.Xor(actual, expected, FillRule.EvenOdd, 4);
        Assert.That(Math.Abs(Clipper.Area(difference)), Is.LessThan(0.01), "Rekeyholing must not retain an old seam");
    }

    [TestCase("ilb_17", 0)]
    [TestCase("ilb_17", 1)]
    [TestCase("ilb_16", 0)]
    [TestCase("ilb_16", 1)]
    [TestCase("ilb_16b", 0)]
    [TestCase("ilb_16b", 1)]
    public void RepeatedBoolean_ReplacesOriginalKeyhole(string project, int legacy)
    {
        Fixture fixture = LoadFixture($"{project}-{legacy}");
        PathsD layerFive = RoundAndBias(fixture, legacy);
        Assert.That(Clipper.Area(layerFive), Is.EqualTo(fixture.ExpectedLayer5Area).Within(0.01));
        PathsD layerSeven = GeoWrangler.customBoolean(0, layerFive, 1, ToPaths(fixture.Clip), 0,
            fixture.Resolution, (double)fixture.BooleanExtension);
        layerSeven = ProcessOverlaps.processOverlaps(layerSeven, layerSeven.Select(_ => false).ToList(),
            (double)fixture.BooleanExtension, fixture.Resolution).geometry;
        foreach (var (points, holeCount) in new[] { (layerFive, 1), (layerSeven, 2) })
        {
            Assert.That(points, Has.Count.EqualTo(1), "Redundant keyholes must not split off a wedge");
            var physical = GeoWrangler.sliverGapRemoval(points, customSizing: 2 * GeoWrangler.keyhole_sizing);
            Assert.That(physical, Has.Count.EqualTo(holeCount + 1), "Preserve the outer and each physical hole");
        }
        Assert.That(Clipper.Area(layerSeven), Is.EqualTo(fixture.ExpectedLayer7Area).Within(0.01));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void Proximity_PreservesInputAndDrawnGeometry(bool drawn)
    {
        PathD outer = Clipper.MakePath(new double[] { 0, 0, 40, 0, 40, 40, 0, 40 });
        PathD hole = Clipper.MakePath(new double[] { 10, 10, 10, 30, 30, 30, 30, 10 });
        PathsD input = GeoWrangler.makeKeyHole(new PathsD { outer }, new PathsD { hole }, false, true);
        foreach (PathD path in input) path.Reverse();
        var original = input.Select(p => p.Select(v => (v.x, v.y)).ToArray()).ToArray();
        var result = Proximity.proximityBias(input, new() { drawn }, 5, 30, 4, 0, 1, 1.03m, 1, false, 0);
        Assert.That(input.Select(p => p.Select(v => (v.x, v.y)).ToArray()).ToArray(), Is.EqualTo(original));
        Assert.That(result.geometry, Is.Not.Empty);
        Assert.That(result.drawn.All(v => v == drawn), Is.True);
        if (drawn)
            Assert.That(Math.Abs(Clipper.Area(result.geometry)), Is.EqualTo(Math.Abs(Clipper.Area(input))).Within(0.001));
    }

    [TestCase(50)]
    [TestCase(30)]
    public void LayerFive_ModernRounding_ProducesCleanKeyhole(int barLength)
    {
        Fixture fixture = LoadFixture($"ilb_4-{barLength}");
        PathsD points = RoundAndBias(fixture, 0);
        Assert.That(Clipper.Area(points), Is.EqualTo(fixture.ExpectedLayer5Area).Within(0.01));
        Assert.That(points, Has.Count.EqualTo(1));
        Assert.That(points.SelectMany(p => p).All(p => double.IsFinite(p.x) && double.IsFinite(p.y)), Is.True);
        PathD boundary = GeoWrangler.stripTerminators(points[0], false);
        for (int i = 0; i < boundary.Count; i++)
        {
            PointD current = boundary[i];
            PointD previous = boundary[(i + boundary.Count - 1) % boundary.Count];
            PointD next = boundary[(i + 1) % boundary.Count];
            double ax = previous.x - current.x, ay = previous.y - current.y;
            double bx = next.x - current.x, by = next.y - current.y;
            double length = Math.Sqrt((ax * ax + ay * ay) * (bx * bx + by * by));
            if (length <= 1e-12) continue;
            double angle = Math.Acos(Math.Clamp((ax * bx + ay * by) / length, -1, 1)) * 180 / Math.PI;
            Assert.That(angle, Is.GreaterThan(45), $"Spike at {current.x},{current.y}");
        }
        if (barLength == 30)
        {
            var contours = GeoWrangler.sliverGapRemoval(points, customSizing: 2 * GeoWrangler.keyhole_sizing);
            Assert.That(contours, Has.Count.EqualTo(2), "The enclosed hole must survive proximity bias");
        }
    }

    private static Fixture LoadFixture(string name) =>
        JsonSerializer.Deserialize<Fixture>(File.ReadAllText(Path.Combine(
            TestContext.CurrentContext.TestDirectory, "TestData", "LayoutRegression", name + ".json")))!;

    private static PathsD ToPaths(double[][][] paths) =>
        new(paths.Select(path => new PathD(path.Select(point => new PointD(point[0], point[1])))));

    // Exercise the library boundary captured from the ILB projects: round the
    // Boolean geometry, rotate, apply proximity bias, then reconcile overlaps.
    private static PathsD RoundAndBias(Fixture fixture, int legacy)
    {
        ShapeSettings settings = new();
        settings.setInt(ShapeSettings.properties_i.enabled, 1);
        settings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        settings.setInt(ShapeSettings.properties_i.legacyRounding, legacy);
        settings.setDecimal(ShapeSettings.properties_decimal.iCR, fixture.InnerRadius);
        settings.setDecimal(ShapeSettings.properties_decimal.oCR, fixture.OuterRadius);
        settings.setDecimal(ShapeSettings.properties_decimal.rot, fixture.Rotation);
        PathsD source = ToPaths(fixture.Source);
        RectD bounds = Clipper.GetBounds(source);
        PointD pivot = new((bounds.left + bounds.right) / 2, (bounds.top + bounds.bottom) / 2);
        Fragmenter fragmenter = new(fixture.Resolution);
        PathsD rounded = [];
        foreach (PathD path in source)
        {
            ShapeLibrary shape = new(Enumerable.Range(0, 12).ToArray(), settings);
            shape.setShape((int)ShapeLibrary.shapeNames_all.GEOCORE, GeoWrangler.clockwiseAndReorderXY(path));
            PathD contour = shape.processCorners(true, false, fixture.CornerSegments,
                fixture.OptimizeCorners, fixture.Resolution);
            contour = ShapeLibrary.rotateShape(contour, settings, 0, 0, pivot).output;
            rounded.Add(fragmenter.fragmentPath(GeoWrangler.close(contour)));
        }
        PathsD sampled = fragmenter.fragmentPaths(GeoWrangler.stripCollinear(rounded));
        PathsD biased = Proximity.proximityBias(sampled, sampled.Select(_ => false).ToList(),
            fixture.Bias, fixture.Distance, fixture.Rays, 0, 1, fixture.Extension,
            fixture.Resolution, false, 0).geometry;
        biased = fragmenter.fragmentPaths(GeoWrangler.stripCollinear(biased));
        return ProcessOverlaps.processOverlaps(biased, biased.Select(_ => false).ToList(),
            (double)fixture.Extension, fixture.Resolution).geometry;
    }

    private sealed class Fixture
    {
        public double[][][] Source { get; set; } = [];
        public double[][][] Clip { get; set; } = [];
        public int CornerSegments { get; set; }
        public int OptimizeCorners { get; set; }
        public double Resolution { get; set; }
        public decimal InnerRadius { get; set; }
        public decimal OuterRadius { get; set; }
        public decimal Rotation { get; set; }
        public double ExpectedLayer5Area { get; set; }
        public double ExpectedLayer7Area { get; set; }
        public decimal Bias { get; set; }
        public decimal Distance { get; set; }
        public int Rays { get; set; }
        public decimal Extension { get; set; }
        public decimal BooleanExtension { get; set; }
    }
}
