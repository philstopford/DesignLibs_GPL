using System;
using System.IO;
using System.Linq;
using geoCoreLib;
using NUnit.Framework;

namespace UnitTests;

public class LayoutImportTests
{
    private static GeoCore Import(string relativePath, GeoCore.fileType format)
    {
        GeoCore geo = new();
        geo.updateGeoCore(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "LayoutRegression", relativePath), format);
        Assert.That(geo.isValid(), Is.True, string.Join(Environment.NewLine, geo.error_msgs));
        return geo;
    }

    [TestCase("doetest.gds", GeoCore.fileType.gds)]
    [TestCase("doetest.oas", GeoCore.fileType.oasis)]
    public void DoeLayout_ImportsInNanometers(string relativePath, GeoCore.fileType format)
    {
        GeoCore geo = Import(relativePath, format);
        Assert.That(geo.getStructureList(), Is.Not.Empty);
        geo.updateGeometry(geo.activeStructure, 0);
        var points = geo.points(true);
        Assert.That(points, Is.Not.Empty);
        double extent = points.SelectMany(p => p).Max(p => Math.Max(Math.Abs(p.x), Math.Abs(p.y)));
        Assert.That(extent, Is.EqualTo(155).Within(1e-9), "DOE geometry should be measured in nanometers");
    }

    [Test]
    public void DoeLayout_RepeatedOasisGeometry_MatchesGds()
    {
        GeoCore gds = Import("doetest.gds", GeoCore.fileType.gds);
        GeoCore oasis = Import("doetest.oas", GeoCore.fileType.oasis);
        Assert.That(oasis.getStructureList(), Is.EquivalentTo(gds.getStructureList()));
        foreach (string structure in gds.getStructureList())
        {
            int gdsIndex = gds.getStructureList().IndexOf(structure);
            int oasisIndex = oasis.getStructureList().IndexOf(structure);
            gds.updateGeometry(gdsIndex, 0);
            oasis.updateGeometry(oasisIndex, 0);
            Assert.That(oasis.getActiveStructureLDList(), Is.EqualTo(gds.getActiveStructureLDList()));
            for (int layer = 0; layer < gds.getActiveStructureLDList().Count; layer++)
            {
                gds.updateGeometry(gdsIndex, layer);
                oasis.updateGeometry(oasisIndex, layer);
                // Compare vertices independently of polygon ordering and starting vertex.
                var expected = gds.points(true).SelectMany(p => p.Select(v => (x: Math.Round(v.x, 9), y: Math.Round(v.y, 9))).Distinct()).OrderBy(p => p.x).ThenBy(p => p.y).ToArray();
                var actual = oasis.points(true).SelectMany(p => p.Select(v => (x: Math.Round(v.x, 9), y: Math.Round(v.y, 9))).Distinct()).OrderBy(p => p.x).ThenBy(p => p.y).ToArray();
                Assert.That(actual, Is.EqualTo(expected), $"{structure}, layer {layer}");
            }
        }
    }
}
