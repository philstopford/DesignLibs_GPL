using Clipper2Lib;
using geoWrangler;
using shapeEngine;

namespace UnitTests;
using System;
using System.Reflection;

// Paste this into your test project (same assembly as ShapeLibrary) and call InvokeLegacyFromTest.
// Requires the project's PathD, ShapeLibrary, GeoWrangler and Clipper types to be in scope.

public static class LegacyReflectionHelper
{
    /// <summary>
    /// Invokes the private legacy_processCorners_actual method on an existing ShapeLibrary instance.
    /// Returns the produced PathD and its signed area (as returned by Clipper.Area).
    /// </summary>
    public static (PathD path, double area) InvokeLegacy(
        ShapeLibrary shapeLib,
        bool previewMode,
        bool cornerCheck,
        int cornerSegments,
        int optimizeCorners,
        double resolution,
        bool iCPA = false,
        bool oCPA = false,
        double iCV = 0.0,
        double iCVariation_scalar = 0.0,
        double oCV = 0.0,
        double oCVariation_scalar = 0.0)
    {
        if (shapeLib == null) throw new ArgumentNullException(nameof(shapeLib));

        var method = typeof(ShapeLibrary).GetMethod(
            "legacy_processCorners_actual",
            BindingFlags.NonPublic | BindingFlags.Instance);

        if (method == null)
            throw new InvalidOperationException("legacy_processCorners_actual method not found on ShapeLibrary");

        object[] args = new object[]
        {
            previewMode,
            cornerCheck,
            cornerSegments,
            optimizeCorners,
            resolution,
            iCPA,
            oCPA,
            iCV,
            iCVariation_scalar,
            oCV,
            oCVariation_scalar
        };

        // If the private method throws, TargetInvocationException will be raised — let it bubble so the test can see the cause.
        var result = method.Invoke(shapeLib, args) as PathD ?? new PathD();

        // Optionally remove duplicate points before area measurement to match other tests:
        var clean = GeoWrangler.removeDuplicates(result);

        double area = Clipper.Area(clean);

        return (clean, area);
    }
}