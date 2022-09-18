using Clipper2Lib;
using geoWrangler;
#pragma warning disable CS8618

namespace geoAnalysis;

using Paths = Paths64;

public class AreaHandler
{
    public enum areaCalcModes { all, perpoly }
    public double area { get; private set; }
    public Paths listOfOutputPoints { get; private set; }

    private void ZFillCallback(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt.Z = -1; // Tag our intersection points.
    }

    public AreaHandler(Paths aPaths, Paths bPaths, bool maySimplify, bool perPoly, double scaleFactorForOperation)
    {
        areaHandlerLogic(aPaths, bPaths, scaleFactorForOperation, maySimplify, perPoly);
    }

    private void areaHandlerLogic(Paths aPaths, Paths bPaths, double scaleFactorForOperation, bool maySimplify, bool perPoly)
    {
        Paths tmpPaths = new();
        listOfOutputPoints = new Paths();

        // callsite may not want simplified geometry.
        Clipper64 c = new() {ZCallback = ZFillCallback, PreserveCollinear = !maySimplify};

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
        area = tmpVal / (scaleFactorForOperation * scaleFactorForOperation);
    }
}