using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static bool isClockwise(Path64 points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(Path64 points)
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
                deltaX = points[0].X - points[pt].X;
                deltaY = points[0].Y + points[pt].Y;
            }
            else
            {
                deltaX = points[pt + 1].X - points[pt].X;
                deltaY = points[pt + 1].Y + points[pt].Y;
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