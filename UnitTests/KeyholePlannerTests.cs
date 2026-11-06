using System;
using System.Linq;
using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;

namespace UnitTests;

public class KeyholePlannerTests
{
    private static PathD Box(double x, double y, double w, double h, bool hole = false)
    {
        PathD path = Clipper.MakePath(new[] { x, y, x + w, y, x + w, y + h, x, y + h });
        return hole ? Clipper.ReversePath(path) : path;
    }

    [TestCase(false, 2.5)]
    [TestCase(true, 2.5)]
    [TestCase(false, 10)]
    [TestCase(true, 10)]
    public void Chain_UsesHoleToHoleConnection_AndIsIndependentOfInputOrder(bool permute, double sizing)
    {
        PathD outer = Box(0, 0, 100, 100);
        PathD first = Box(30, 40, 10, 10, true), second = Box(45, 40, 10, 10, true);
        PathsD holes = permute ? [second, first] : [first, second];
        var snapshot = holes.Select(p => p.ToArray()).ToArray();
        PathsD result = GeoWrangler.makeKeyHole(new PathsD { outer }, holes, false, true, customSizing: sizing);
        Assert.That(holes.Select(p => p.ToArray()).ToArray(), Is.EqualTo(snapshot));
        Assert.That(result, Has.Count.EqualTo(1));
        // MST: outer -> first = 30, first -> second = 5. Separate
        // outer connections would total 70. Corridor width is 0.1;
        // three exposed square-cap corners add 3 * 0.0025.
        Assert.That(9800 - Math.Abs(Clipper.Area(result)), Is.EqualTo(35 * (0.04 * sizing) + 3 * Math.Pow(0.02 * sizing, 2)).Within(0.0001));
        PathsD reference = GeoWrangler.makeKeyHole(new PathsD { outer }, new PathsD { first, second }, false, true, customSizing: sizing);
        Assert.That(Math.Abs(Clipper.Area(Clipper.Xor(result, reference, FillRule.EvenOdd, 4))), Is.LessThan(0.0001));
        Assert.That(GeoWrangler.sliverGapRemoval(result, customSizing: sizing), Has.Count.EqualTo(3));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void NestedIslands_AndSeparateComponents_PreserveAllMaterial(bool separateInputs)
    {
        PathsD input = [Box(0, 0, 100, 100), Box(10, 10, 80, 80, true),
            Box(20, 20, 60, 60), Box(30, 30, 40, 40, true),
            Box(200, 0, 100, 100), Box(220, 20, 20, 20, true)];
        var snapshot = input.Select(p => p.ToArray()).ToArray();
        PathsD result = separateInputs
            ? GeoWrangler.makeKeyHole(new PathsD(input.Where(Clipper.IsPositive)),
                new PathsD(input.Where(p => !Clipper.IsPositive(p))), false, true)
            : GeoWrangler.makeKeyHole(input, false, true);
        Assert.That(result, Has.Count.EqualTo(3));
        Assert.That(input.Select(p => p.ToArray()).ToArray(), Is.EqualTo(snapshot));
        PathsD physical = GeoWrangler.sliverGapRemoval(result);
        Assert.That(physical, Has.Count.EqualTo(6));
        Assert.That(Math.Abs(Clipper.Area(physical)), Is.EqualTo(15200).Within(0.05));
    }

    [Test]
    public void ConcaveBoundary_RejectsConnectionsThroughExterior()
    {
        PathD outer = Clipper.MakePath(new double[] { 0,0, 100,0, 100,100, 60,100, 60,20, 40,20, 40,100, 0,100 });
        PathsD holes = [Box(20, 70, 10, 10, true), Box(70, 70, 10, 10, true)];
        PathsD result = GeoWrangler.makeKeyHole(new PathsD { outer }, holes, false, true);
        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(Math.Abs(Clipper.Area(result)), Is.EqualTo(8198).Within(0.01));
        Assert.That(GeoWrangler.sliverGapRemoval(result), Has.Count.EqualTo(3));
    }

    [Test]
    public void OverlappingCutters_AreUnionedWithoutLeavingAnIsland()
    {
        PathsD result = GeoWrangler.makeKeyHole(new PathsD { Box(0, 0, 100, 100) },
            new PathsD { Box(30, 30, 20, 20, true), Box(40, 40, 20, 20, true) }, false, true);
        Assert.That(result, Has.Count.EqualTo(1));
        PathsD physical = GeoWrangler.sliverGapRemoval(result);
        Assert.That(physical, Has.Count.EqualTo(2));
        Assert.That(Math.Abs(Clipper.Area(physical)), Is.EqualTo(9300).Within(0.01));
        Assert.That(GeoWrangler.makeKeyHole(new PathsD(), new PathsD { Box(0, 0, 1, 1, true) }, false, true), Is.Empty);
        Assert.That(GeoWrangler.makeKeyHole(new PathsD { Box(0, 0, 10, 10) },
            new PathsD { Box(-1, -1, 12, 12, true) }, false, true), Is.Empty);
    }

    [Test]
    public void ImpossibleCorridorWidth_FailsWithoutReturningPartialGeometry()
    {
        Assert.Throws<InvalidOperationException>(() => GeoWrangler.makeKeyHole(
            new PathsD { Box(0, 0, 10, 10) }, new PathsD { Box(4, 4, 2, 2, true) },
            false, true, customSizing: 1000));
    }

    [Test]
    public void ReversedWinding_AndRepeatedVertices_ChooseTheSameTree()
    {
        PathsD source = [Box(0, 0, 100, 100), Box(30, 40, 10, 10, true), Box(45, 40, 10, 10, true)];
        PathsD reference = GeoWrangler.makeKeyHole(source, false, true);
        foreach (PathD path in source)
        {
            path.Reverse();
            path.Insert(1, path[0]);
            path.Add(path[0]);
        }
        source.Reverse();
        PathsD result = GeoWrangler.makeKeyHole(source, false, true);
        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(Math.Abs(Clipper.Area(Clipper.Xor(result, reference, FillRule.EvenOdd, 4))), Is.LessThan(0.0001));
    }

    [Test]
    public void EmptyAndDegenerateBoundaries_DoNotRequireConnections()
    {
        Assert.That(GeoWrangler.makeKeyHole(new PathsD(), false, true), Is.Empty);
        PathsD result = GeoWrangler.makeKeyHole(new PathsD { Box(0, 0, 10, 10) },
            new PathsD { Clipper.MakePath(new double[] { 2,2, 3,2, 4,2 }) }, false, true);
        Assert.That(result, Has.Count.EqualTo(1));
        Assert.That(Math.Abs(Clipper.Area(result)), Is.EqualTo(100));
    }
}
