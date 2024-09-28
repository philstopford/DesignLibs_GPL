using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathD getBounds(PathD source)
    {
        return pGetBounds(source);
    }

    private static PathD pGetBounds(PathD source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count <= 0)
        {
            return [new PointD(minX, minY), new PointD(maxX, maxY)];
        }
        try
        {
            minX = source[MinX(source)].x;
            maxX = source[MaxX(source)].x;
            minY = source[MinY(source)].y;
            maxY = source[MaxY(source)].y;
        }
        catch
        {
            // ignored
        }

        return [new PointD(minX, minY), new PointD(maxX, maxY)];
    }
    
    public static Path64 getBounds(Path64 source)
    {
        return pGetBounds(source);
    }

    private static Path64 pGetBounds(Path64 source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count <= 0)
        {
            return [new Point64(minX, minY), new Point64(maxX, maxY)];
        }
        try
        {
            minX = source[MinX(source)].X;
            maxX = source[MaxX(source)].X;
            minY = source[MinY(source)].Y;
            maxY = source[MaxY(source)].Y;
        }
        catch
        {
            // ignored
        }

        return [new Point64(minX, minY), new Point64(maxX, maxY)];
    }
}