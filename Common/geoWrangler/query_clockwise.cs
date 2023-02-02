using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static bool isClockwise(PathD points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(PathD points)
    {
        // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        // Shoelace formula.

        double delta = 0;

        for (int pt = 0; pt < points.Count; pt++)
        {
            double deltaX;
            double deltaY;
            if (pt == points.Count - 1)
            {
                deltaX = points[0].x - points[pt].x;
                deltaY = points[0].y + points[pt].y;
            }
            else
            {
                deltaX = points[pt + 1].x - points[pt].x;
                deltaY = points[pt + 1].y + points[pt].y;
            }

            delta += deltaX * deltaY;
        }

        return delta switch
        {
            > 0 => true,
            _ => false
        };
    }
}