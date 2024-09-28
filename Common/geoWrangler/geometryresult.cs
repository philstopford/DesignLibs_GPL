using System.Collections.Generic;
using Clipper2Lib;

namespace geoWrangler;

public class GeometryResult
{
    public PathsD geometry { get; init; } = [];
    public List<bool> drawn { get; init; } = [];
}
