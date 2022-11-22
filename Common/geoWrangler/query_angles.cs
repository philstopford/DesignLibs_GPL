using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static double[] angles(PathD sourcePoly, bool allowNegative)
    {
        return pAngles(sourcePoly, allowNegative);
    }

    private static double[] pAngles(PathD sourcePoly, bool allowNegative)
    {
        PathD stripped = pStripTerminators(sourcePoly, false);
        int finalIndex = stripped.Count - 1;

        double[] angles = new double[stripped.Count];

        for (int pt = 0; pt <= finalIndex; pt++)
        {
            // Assess angle.
            PointD interSection_A;
            PointD interSection_B;
            PointD interSection_C;
            switch (pt)
            {
                case 0:
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                    break;
                default:
                {
                    if (pt == finalIndex) // last point in the list
                    {
                        interSection_B = stripped[pt - 1];
                        interSection_C = stripped[pt];
                        interSection_A = stripped[0]; // map to the first point
                    }
                    else
                    {
                        interSection_B = stripped[pt - 1];
                        interSection_C = stripped[pt];
                        interSection_A = stripped[pt + 1];
                    }

                    break;
                }
            }

            double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);

            angles[pt] = theta;
        }

        return angles;
    }
}