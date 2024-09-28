using Clipper2Lib;

namespace geoLib;

public class GeoLibArray
{
    public Point64 point { get; set; }
    public Point64 pitch { get; init; }
    public Point64 count { get; init; }
}

public class GeoLibArrayF
{
    public PointD point { get; set; }
    public PointD pitch { get; set; }
    public Point64 count { get; set; }
}