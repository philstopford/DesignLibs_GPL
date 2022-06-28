using System.Collections.Generic;
using geoLib;

namespace geoWrangler;

public class GeometryResult
{
    public List<GeoLibPointF[]> geometry { get; set; }
    public List<bool> drawn { get; set; }

    public GeometryResult()
    {
        geometry = new();
        drawn = new();
    }
}
