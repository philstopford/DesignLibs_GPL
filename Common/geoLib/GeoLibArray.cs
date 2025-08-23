using Clipper2Lib;

namespace geoLib;

/// <summary>
/// Represents an array configuration with integer coordinates for geometric operations.
/// </summary>
public class GeoLibArray
{
    /// <summary>Base point for the array</summary>
    public Point64 point { get; set; }
    
    /// <summary>Spacing between array elements (read-only after initialization)</summary>
    public Point64 pitch { get; init; }
    
    /// <summary>Number of elements in each direction (read-only after initialization)</summary>
    public Point64 count { get; init; }
}

/// <summary>
/// Represents an array configuration with floating-point coordinates for geometric operations.
/// </summary>
public class GeoLibArrayF
{
    /// <summary>Base point for the array (floating-point)</summary>
    public PointD point { get; set; }
    
    /// <summary>Spacing between array elements (floating-point)</summary>
    public PointD pitch { get; set; }
    
    /// <summary>Number of elements in each direction</summary>
    public Point64 count { get; set; }
}