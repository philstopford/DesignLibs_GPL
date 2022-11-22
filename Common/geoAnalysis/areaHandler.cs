using Clipper2Lib;
using geoWrangler;
#pragma warning disable CS8618

namespace geoAnalysis;

public class AreaHandler
{
    public enum areaCalcModes { all, perpoly }
    public double area { get; private set; }
    public PathsD listOfOutputPoints { get; private set; }

    private void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = -1; // Tag our intersection points.
    }

    public AreaHandler(PathsD aPaths, PathsD bPaths, bool maySimplify, bool perPoly)
    {
        areaHandlerLogic(aPaths, bPaths, maySimplify, perPoly);
    }

    private void areaHandlerLogic(PathsD aPaths, PathsD bPaths, bool maySimplify, bool perPoly)
    {
        PathsD tmpPaths = new();
        listOfOutputPoints = new ();

        // callsite may not want simplified geometry.
        ClipperD c = new() {ZCallback = ZFillCallback, PreserveCollinear = !maySimplify};

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
                listOfOutputPoints.Add(tmpPaths[poly]);
            }
            else
            {
                tmpVal += Clipper.Area(tmpPaths[poly]);
                // Append the result output to the resultPoints list.
                listOfOutputPoints.Add(tmpPaths[poly]);
            }
        }
        // Sum the areas by polygon.
        area = tmpVal;
    }
}