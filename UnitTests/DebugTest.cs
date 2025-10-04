using NUnit.Framework;
using Clipper2Lib;
using geoWrangler;
using shapeEngine;
using UnitTests;

public class DebugTest
{
    [Test]
    public static void DebugCustomOrthoOuterRoundingTest()
    {
        ShapeSettings shapeSettings = new ShapeSettings();
        shapeSettings.setInt(ShapeSettings.properties_i.shapeIndex, (int)ShapeLibrary.shapeNames_all.GEOCORE);
        shapeSettings.setInt(ShapeSettings.properties_i.legacyRounding, 1);
        PathD customShape = Clipper.MakePath(new double[]
        {
            0, 0,
            0, 20,
            20, 20,
            20, 0
        });
        shapeSettings.setDecimal(ShapeSettings.properties_decimal.oCR, 10);
        ShapeLibrary shape = new ShapeLibrary(ShapeEngineTests.shapeTable, shapeSettings);
        shape.setShape(shapeSettings.getInt(ShapeSettings.properties_i.shapeIndex), customShape);
        
        PathD out_ = shape.processCorners(false, false, 90, 1, .01);
        PathD clean = GeoWrangler.removeDuplicates(out_);
        double area = Clipper.Area(out_);
        
        var (legacyPath, legacyArea) = LegacyReflectionHelper.InvokeLegacy(
            shape,
            previewMode: false,
            cornerCheck: false,
            cornerSegments: 90,
            optimizeCorners: 1,
            resolution: 0.01,
            iCPA: false,
            oCPA: false,
            iCV: 0, iCVariation_scalar: 0, oCV: 0, oCVariation_scalar: 0
        );
        
        Console.WriteLine($"processCorners output:");
        Console.WriteLine($"  clean.Count = {clean.Count}");
        Console.WriteLine($"  area = {area}");
        Console.WriteLine($"Legacy output:");
        Console.WriteLine($"  legacyPath.Count = {legacyPath.Count}");
        Console.WriteLine($"  legacyArea = {legacyArea}");
        Console.WriteLine($"Expected area: {-(Math.PI * 10 * 10)}");
    }
}
