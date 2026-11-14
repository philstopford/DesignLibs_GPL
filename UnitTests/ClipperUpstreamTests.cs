using System.Text.Json;
using Clipper2Lib;
using NUnit.Framework;

namespace UnitTests;

[TestFixture]
public class ClipperUpstreamTests
{
    // Native reference outputs from PolyClipping/Clipper2 884ba82 (USINGZ, default precision).
    public static IEnumerable<TestCaseData> NativeCases()
    {
        foreach (string file in new[] { "polygons", "lines" })
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Clipper2", file + ".json");
            using JsonDocument doc = JsonDocument.Parse(File.ReadAllText(path));
            foreach (JsonElement entry in doc.RootElement.EnumerateArray())
                yield return new TestCaseData(entry.Clone()).SetName($"Native_{file}_{entry.GetProperty("number").GetInt32()}");
        }
    }

    public static IEnumerable<TestCaseData> NativeOffsets()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Clipper2", "offsets.json");
        using JsonDocument doc = JsonDocument.Parse(File.ReadAllText(path));
        foreach (JsonElement entry in doc.RootElement.EnumerateArray())
            yield return new TestCaseData(entry.Clone()).SetName($"Native_offset_{entry.GetProperty("number").GetInt32()}");
    }

    [TestCaseSource(nameof(NativeOffsets))]
    public void OffsettingMatchesNativeReference(JsonElement entry)
    {
        // C++ orders the JoinType enum differently from the established managed API.
        JoinType[] nativeJoins = [JoinType.Square, JoinType.Bevel, JoinType.Round, JoinType.Miter];
        var co = new ClipperOffset(2, 0.25) { Compatibility = ClipperCompatibility.Upstream };
        co.AddPaths(ReadPaths(entry, "subject"), nativeJoins[entry.GetProperty("joinType").GetInt32()],
            (EndType)entry.GetProperty("endType").GetInt32());
        var result = new Paths64();
        double delta = entry.GetProperty("delta").GetDouble();
        co.Execute(delta, result);
        Assert.That(Canonical(result), Is.EqualTo(Canonical(ReadPaths(entry, "result"))));
        var tree = new PolyTree64(); co.Execute(delta, tree);
        Assert.That(Canonical(Clipper.PolyTreeToPaths64(tree)), Is.EqualTo(Canonical(result)));
    }

    private static Paths64 ReadPaths(JsonElement entry, string name) => new(entry.GetProperty(name).EnumerateArray()
        .Select(p => new Path64(p.EnumerateArray().Select(pt => new Point64(pt[0].GetInt64(), pt[1].GetInt64())))));

    // Ignore path list order and the arbitrary starting vertex, but preserve winding and all vertices.
    private static string[] Canonical(Paths64 paths, bool open = false) => paths.Select(path =>
    {
        if (path.Count == 0) return "";
        string Vertex(int i) => $"{path[i].X},{path[i].Y}";
        if (open) return string.Join(";", Enumerable.Range(0, path.Count).Select(Vertex));
        return Enumerable.Range(0, path.Count).Select(start => string.Join(";",
            Enumerable.Range(0, path.Count).Select(i => Vertex((start + i) % path.Count))))
            .Min(StringComparer.Ordinal)!;
    }).Order(StringComparer.Ordinal).ToArray();

    [TestCaseSource(nameof(NativeCases))]
    public void ClippingMatchesNativeReference(JsonElement entry)
    {
        Paths64 subject = ReadPaths(entry, "subject"), open = ReadPaths(entry, "open"), clip = ReadPaths(entry, "clip");
        Paths64 expected = ReadPaths(entry, "closedResult"), expectedOpen = ReadPaths(entry, "openResult");
        var c = new Clipper64 { Compatibility = ClipperCompatibility.Upstream };
        c.AddSubject(subject); c.AddOpenSubject(open); c.AddClip(clip);
        var result = new Paths64(); var resultOpen = new Paths64();
        var ct = (ClipType)entry.GetProperty("clipType").GetInt32();
        var fr = (FillRule)entry.GetProperty("fillRule").GetInt32();
        Assert.That(c.Execute(ct, fr, result, resultOpen), Is.True);
        Assert.That(Canonical(result), Is.EqualTo(Canonical(expected)));
        Assert.That(Canonical(resultOpen, true), Is.EqualTo(Canonical(expectedOpen, true)));
        // Force both nesting strategies, then repeat on the same instance to exercise pooled state.
        foreach (int threshold in new[] { 0, int.MaxValue, 64 })
        {
            c.PipIndexThreshold = threshold;
            var tree = new PolyTree64();
            Assert.That(c.Execute(ct, fr, tree, resultOpen), Is.True);
            Assert.That(Canonical(Clipper.PolyTreeToPaths64(tree)), Is.EqualTo(Canonical(expected)));
            Assert.That(Canonical(resultOpen, true), Is.EqualTo(Canonical(expectedOpen, true)));
        }
    }

    private static Path64 P(params long[] coords) => Clipper.MakePath(coords);
    private static Path64 Square(long x, long y, long size, bool clockwise = false)
    {
        var p = P(x, y, x + size, y, x + size, y + size, x, y + size);
        if (clockwise) p.Reverse();
        return p;
    }
    private static Paths64 Offset(Paths64 paths, double delta, JoinType jt = JoinType.Miter, EndType et = EndType.Polygon)
    {
        var co = new ClipperOffset { PreserveInputOrientation = false, Compatibility = ClipperCompatibility.Upstream };
        co.AddPaths(paths, jt, et);
        var result = new Paths64(); co.Execute(delta, result); return result;
    }

    [TestCase(false)]
    [TestCase(true)]
    public void OffsetRetainsLegacyWindingAndOffersUpstreamWinding(bool reverse)
    {
        Path64 square = Square(0, 0, 100, true);
        var co = new ClipperOffset { ReverseSolution = reverse };
        co.AddPath(square, JoinType.Miter, EndType.Polygon);
        var result = new Paths64(); co.Execute(10, result);
        Assert.That(Clipper.Area(result), Is.EqualTo(reverse ? 14400 : -14400));
        Assert.That(Clipper.Area(Clipper.InflatePaths(new Paths64 { square }, 10, JoinType.Miter, EndType.Polygon)), Is.EqualTo(-14400));
        co.PreserveInputOrientation = false;
        co.Compatibility = ClipperCompatibility.Upstream;
        co.Execute(10, result);
        Assert.That(Clipper.Area(result), Is.EqualTo(reverse ? -14400 : 14400));
    }

    [TestCase(0.0)]
    [TestCase(10.0)]
    [TestCase(-10.0)]
    public void OffsetNestedGroupsAreIndependentOfInsertionOrder(double delta)
    {
        Path64[] chain = [Square(0, 0, 1000), Square(200, 200, 500, true), Square(300, 300, 100)];
        var expected = Offset(new Paths64(chain), delta, JoinType.Round);
        foreach (int[] order in new[] { new[] { 0, 1, 2 }, [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0] })
        {
            var co = new ClipperOffset { PreserveInputOrientation = false, Compatibility = ClipperCompatibility.Upstream };
            foreach (int i in order) co.AddPath(chain[i], JoinType.Round, EndType.Polygon);
            var result = new Paths64(); co.Execute(delta, result);
            Assert.That(Canonical(result), Is.EqualTo(Canonical(expected)));
        }
    }

    [Test]
    public void OffsetJoinedShortPathDoesNotChangeLaterCaps()
    {
        var two = P(0, 0, 100, 0); var tri = P(0, 1000, 500, 1000, 500, 1500);
        var result = Offset(new Paths64 { two, tri }, 20, JoinType.Miter, EndType.Joined);
        Assert.That(Clipper.Area(result), Is.EqualTo(73268.5));
        Assert.That(Canonical(result), Is.EqualTo(Canonical(Offset(new Paths64 { tri, two }, 20, JoinType.Miter, EndType.Joined))));
    }

    [TestCase(10.0)]
    [TestCase(-10.0)]
    public void OffsetDropsStrayHolesAndUsesExactArea(double delta)
    {
        Path64 outer = Square(0, 0, 100), stray = Square(500, 500, 50, true);
        Assert.That(Canonical(Offset(new Paths64 { outer, stray }, delta)), Is.EqualTo(Canonical(Offset(new Paths64 { outer }, delta))));
        var thinOuter = P(-310414360,-251832566, 13677104,72258904, 337768562,396350368, 13677098,72258898);
        var thinHole = P(13677099,72258899, 229738075,288319879, 13677103,72258903, -202383873,-143802077);
        var co = new ClipperOffset { PreserveInputOrientation = false, Compatibility = ClipperCompatibility.Upstream };
        co.AddPath(Square(-1000000, -1000000, 100000), JoinType.Miter, EndType.Polygon);
        co.AddPaths(new Paths64 { thinOuter, thinHole }, JoinType.Miter, EndType.Polygon);
        var result = new Paths64(); co.Execute(delta, result);
        var all = Offset(new Paths64 { Square(-1000000, -1000000, 100000), thinOuter, thinHole }, delta);
        Assert.That(Canonical(result), Is.EqualTo(Canonical(all)));
    }

    [Test]
    public void OffsetCallbackGetsOwnNormalsAndArcSteps()
    {
        var point = P(0, 0); var line = P(0, 1000, 100, 1000);
        Paths64 Run(Paths64 paths)
        {
            var co = new ClipperOffset(); co.AddPaths(paths, JoinType.Round, EndType.Joined);
            var result = new Paths64();
            co.Execute((path, normals, _, _) => { Assert.That(normals.Count, Is.EqualTo(path.Count)); return 20; }, result);
            return result;
        }
        Assert.That(Canonical(Run(new Paths64 { point, line })), Is.EqualTo(Canonical(Run(new Paths64 { line, point }))));
    }

    [Test]
    public void OffsetCanAlternateTreeAndPathOutputs()
    {
        var co = new ClipperOffset(); co.AddPath(Square(0, 0, 100), JoinType.Miter, EndType.Polygon);
        var tree = new PolyTree64(); var paths = new Paths64();
        co.Execute(10, tree); co.Execute(20, paths);
        Assert.That(Clipper.Area(paths), Is.EqualTo(19600));
        Assert.That(tree.Area(), Is.EqualTo(14400));
        co.Execute(5, tree); co.Execute(10, paths);
        Assert.That(Clipper.Area(paths), Is.EqualTo(14400));
    }

    [Test]
    public void EmptyAndDegenerateOffsetsDoNotChangeLaterGroups()
    {
        var square = Square(0, 0, 100);
        var co = new ClipperOffset();
        co.AddPaths(new Paths64 { new(), P(0, 500, 100, 500) }, JoinType.Miter, EndType.Polygon);
        co.AddPath(square, JoinType.Miter, EndType.Polygon);
        var result = new Paths64(); co.Execute(-10, result);
        Assert.That(Clipper.Area(result), Is.EqualTo(6400));
        Assert.That(Offset(new Paths64 { P(500, 500) }, -10), Is.Empty);
        Assert.That(Offset(new Paths64 { P(0, 0, 100, 0) }, -10), Is.Empty);
        Assert.That(Offset(new Paths64 { P(200, 0, 300, 0, 300, 100) }, 0.3, et: EndType.Joined), Is.Empty);
    }

    [Test]
    public void StripeIndexMatchesBruteForceAndCanBeRebuilt()
    {
        var random = new Random(410);
        var index = new StripeIndex();
        foreach (Path64 path in new[] { Clipper.Ellipse(new Point64(0, 0), 1000, 800, 200),
                     P(0,0, 100,0, 100,40, 40,40, 40,100, 0,100), P(0,0, 10,0, 20,0), new Path64() })
        {
            index.BuildFromPath(path);
            if (index.Empty) continue;
            foreach (Point64 pt in path.Concat(Enumerable.Range(0, 2000).Select(_ => new Point64(random.Next(-1100, 1100), random.Next(-900, 900)))))
                Assert.That(index.PointInPolygon(pt), Is.EqualTo(Clipper.PointInPolygon(pt, path)), $"Point {pt}");
        }
        Assert.That(index.Empty, Is.True);
    }

    [Test]
    public void IntersectionUsesEvenRoundingAndIgnoresSegmentDirection()
    {
        var a = new Point64(-207,-383); var b = new Point64(350,194);
        var c = new Point64(347,193); var d = new Point64(-206,-384);
        Assert.That(InternalClipper.GetLineIntersectPt(a, b, c, d, out var pt), Is.True);
        Assert.That(pt, Is.EqualTo(new Point64(66, -100)));
        InternalClipper.GetLineIntersectPt(b, a, d, c, out var reversed);
        Assert.That(reversed, Is.EqualTo(pt));
        var modern = new Clipper64 { Compatibility = ClipperCompatibility.Upstream };
        modern.AddSubject(P(0,10, 5,2, 7,1)); modern.AddClip(P(3,6, 6,1, 11,11));
        var diff = new Paths64(); modern.Execute(ClipType.Difference, FillRule.NonZero, diff);
        Assert.That(diff, Has.Count.EqualTo(1));
        Assert.That(Clipper.Area(diff), Is.EqualTo(3.27).Within(1));
    }

    [TestCase(JoinType.Square)]
    [TestCase(JoinType.Round)]
    [TestCase(JoinType.Bevel)]
    [TestCase(JoinType.Miter)]
    public void OffsetConstructedVerticesKeepZ(JoinType join)
    {
        var path = new Path64 { new(0,0,42), new(100,0,42), new(100,100,42) };
        var co = new ClipperOffset { Compatibility = ClipperCompatibility.Upstream };
        co.AddPath(path, join, EndType.Square);
        var result = new Paths64(); co.Execute(10, result);
        Assert.That(result, Is.Not.Empty);
        Assert.That(result.SelectMany(p => p).All(pt => pt.Z == 42), Is.True);
    }

    [Test]
    public void ScanbeamClampDoesNotLeaveDistantIntersection()
    {
        var co = new Clipper64 { Compatibility = ClipperCompatibility.Upstream };
        co.AddSubject(P(-69425,73, -33569,-35, -24000,-87));
        co.AddClip(P(-3862,183, -33569,-35, -28785,-61, -26393,-75, -24000,-87));
        var result = new Paths64(); co.Execute(ClipType.Union, FillRule.EvenOdd, result);
        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(result[0], Has.Count.EqualTo(4));
        Assert.That(Clipper.Area(result), Is.EqualTo(2230933));
    }

    [Test]
    public void ExactProductComparisonHandlesSignedExtremes()
    {
        Assert.That(InternalClipper.CompareProducts(long.MinValue, long.MinValue, long.MaxValue, long.MaxValue), Is.EqualTo(1));
        Assert.That(InternalClipper.CompareProducts(long.MinValue, 1, long.MaxValue, -1), Is.EqualTo(-1));
        Assert.That(InternalClipper.CompareProducts(-1, -1, 1, 1), Is.Zero);
        Assert.That(InternalClipper.CompareProducts(0, long.MinValue, 0, long.MaxValue), Is.Zero);
    }

    [Test]
    public void ThrowingZCallbackLeavesClipperReusable()
    {
        var c = new Clipper64 { Compatibility = ClipperCompatibility.Upstream };
        c.AddSubject(new Paths64 { Square(0,0,100), Square(50,50,100) });
        c.AddClip(P(20,-10, 80,-10, 80,200, 20,200));
        var expected = new Paths64(); c.Execute(ClipType.Xor, FillRule.NonZero, expected);
        int calls = 0;
        c.ZCallback = (Point64 _, Point64 _, Point64 _, Point64 _, ref Point64 _) =>
        { if (++calls == 3) throw new InvalidOperationException("callback"); };
        var result = new Paths64();
        // Preserve the managed port's boolean failure contract.
        Assert.That(c.Execute(ClipType.Xor, FillRule.NonZero, result), Is.False);
        c.ZCallback = null;
        Assert.That(c.Execute(ClipType.Xor, FillRule.NonZero, result), Is.True);
        Assert.That(Canonical(result), Is.EqualTo(Canonical(expected)));
    }
}
