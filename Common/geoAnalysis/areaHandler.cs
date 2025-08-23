using Clipper2Lib;
using geoWrangler;
#pragma warning disable CS8618

namespace geoAnalysis;

/// <summary>
/// Calculates intersection areas between two sets of polygonal paths.
/// This class performs boolean operations to determine overlapping regions
/// and computes their total or per-polygon areas.
/// </summary>
public class AreaHandler
{
    /// <summary>
    /// Defines calculation modes for area analysis.
    /// </summary>
    public enum areaCalcModes 
    { 
        /// <summary>Calculate total area of all intersections</summary>
        all, 
        /// <summary>Calculate area of individual polygons separately</summary>
        perpoly 
    }
    
    /// <summary>
    /// Gets the calculated intersection area. In 'all' mode, this is the total area.
    /// In 'perpoly' mode, this is the area of the smallest intersecting polygon.
    /// </summary>
    public double area { get; private set; }
    
    /// <summary>
    /// Gets the list of intersection polygons resulting from the boolean operation.
    /// These paths represent the geometric areas where input polygons overlap.
    /// </summary>
    public PathsD listOfOutputPoints { get; private set; }

    /// <summary>
    /// Callback function for Z-fill operations during polygon clipping.
    /// Tags intersection points with a specific Z value for identification.
    /// </summary>
    /// <param name="bot1">Bottom point of first edge</param>
    /// <param name="top1">Top point of first edge</param>
    /// <param name="bot2">Bottom point of second edge</param>
    /// <param name="top2">Top point of second edge</param>
    /// <param name="pt">Reference to the intersection point to be tagged</param>
    private static void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = -1; // Tag our intersection points.
    }

    /// <summary>
    /// Initializes a new instance of the AreaHandler class and performs area calculation.
    /// </summary>
    /// <param name="aPaths">First set of polygonal paths (subject polygons)</param>
    /// <param name="bPaths">Second set of polygonal paths (clip polygons)</param>
    /// <param name="maySimplify">If true, allows simplification of collinear points during processing</param>
    /// <param name="perPoly">If true, calculates per-polygon areas; if false, calculates total area</param>
    public AreaHandler(PathsD aPaths, PathsD bPaths, bool maySimplify, bool perPoly)
    {
        areaHandlerLogic(aPaths, bPaths, maySimplify, perPoly);
    }

    /// <summary>
    /// Core logic for performing area calculations between polygon sets.
    /// Uses boolean intersection operations to find overlapping regions.
    /// </summary>
    /// <param name="aPaths">Subject polygons for intersection</param>
    /// <param name="bPaths">Clip polygons for intersection</param>
    /// <param name="maySimplify">Whether to allow polygon simplification</param>
    /// <param name="perPoly">Whether to calculate per-polygon or total area</param>
    private void areaHandlerLogic(PathsD aPaths, PathsD bPaths, bool maySimplify, bool perPoly)
    {
        PathsD tmpPaths = [];
        listOfOutputPoints = [];

        // callsite may not want simplified geometry.
        ClipperD c = new(Constants.roundingDecimalPrecision) {ZCallback = ZFillCallback, PreserveCollinear = !maySimplify};

        c.AddSubject(aPaths);

        c.AddClip(bPaths);

        // Boolean AND of the two levels for the area operation.
        try
        {
            c.Execute(ClipType.Intersection, FillRule.EvenOdd, tmpPaths); //, firstLayerFillType, secondLayerFillType);
            tmpPaths = GeoWrangler.reOrderXY(tmpPaths);
        }
        catch (Exception)
        {
            // Will handle downstream.
        }
        
        double tmpVal = 0.0;
        if (perPoly)
        {
            tmpVal = -1.0f;
        }

        int polyCount = tmpPaths.Count;
        for (int poly = 0; poly < polyCount; poly++)
        {
            if (perPoly)
            {
                double tmpVal2 = Clipper.Area(tmpPaths[poly]);
                if (!(tmpVal <= -0.0001f) && !(tmpVal2 < tmpVal))
                {
                    continue;
                }

                tmpVal = tmpVal2;
                listOfOutputPoints.Clear();
            }
            else
            {
                tmpVal += Clipper.Area(tmpPaths[poly]);
                // Append the result output to the resultPoints list.
            }

            listOfOutputPoints.Add(tmpPaths[poly]);
        }
        // Sum the areas by polygon.
        area = tmpVal;
    }
}