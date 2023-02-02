using System.Collections.Generic;
using Clipper2Lib;

namespace geoWrangler;

public class GeometryResult
{
    public PathsD geometry { get; set; }
    public List<bool> drawn { get; set; }

    public GeometryResult()
    {
        geometry = new ();
        drawn = new List<bool>();
    }
}
