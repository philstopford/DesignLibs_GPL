using Clipper2Lib;

namespace geoLib;

/// <summary>Direction type for geometry elements</summary>
public enum typeDirection { left1, right1, up1, down1, tilt1 }

/// <summary>Vertex type classification</summary>
public enum typeVertex { corner, center }

/// <summary>Rounding type for geometric elements</summary>
public enum typeRound { inner, exter }

/// <summary>
/// Represents a vertex in 2D space with additional geometric properties.
/// Used for defining geometric shapes with direction, type, and bias tracking.
/// </summary>
public class MyVertex
{
    private PointD internal_pt;

    /// <summary>Gets or sets the X coordinate of the vertex</summary>
    public double X
    {
        get => internal_pt.x;
        set => internal_pt.x = value;
    }
    
    /// <summary>Gets or sets the Y coordinate of the vertex</summary>
    public double Y
    {
        get => internal_pt.y;
        set => internal_pt.y = value;
    }

    /// <summary>Direction associated with this vertex</summary>
    public typeDirection direction { get; set; }
    
    /// <summary>Indicates if this vertex is attached to a vertical edge</summary>
    public bool vertical { get; set; }
    
    /// <summary>Indicates if this is an inner vertex</summary>
    public bool inner { get; set; }
    
    /// <summary>Defines whether the point is a vertex or a center point on an edge</summary>
    public typeVertex type { get; set; }

    /// <summary>Tracks whether X bias has been applied to this vertex</summary>
    public bool xBiasApplied { get; set; }
    
    /// <summary>Tracks whether Y bias has been applied to this vertex</summary>
    public bool yBiasApplied { get; set; }

    /// <summary>
    /// Creates a new vertex with specified coordinates and properties.
    /// </summary>
    /// <param name="X">X coordinate</param>
    /// <param name="Y">Y coordinate</param>
    /// <param name="direction">Direction type</param>
    /// <param name="vertical">True if attached to vertical edge</param>
    /// <param name="inner">True if inner vertex</param>
    /// <param name="type">Vertex type (corner or center)</param>
    public MyVertex(double X, double Y, typeDirection direction, bool vertical, bool inner, typeVertex type)
    {
        pMyVertex(X, Y, direction, vertical, inner, type);
    }

    /// <summary>
    /// Creates a copy of an existing vertex.
    /// </summary>
    /// <param name="source">Source vertex to copy from</param>
    public MyVertex(MyVertex source)
    {
        pMyVertex(source.X, source.Y, source.direction, source.vertical, source.inner, source.type);
    }

    /// <summary>
    /// Internal initialization method for vertex properties.
    /// </summary>
    private void pMyVertex(double X_, double Y_, typeDirection direction_, bool vertical_, bool inner_, typeVertex type_)
    {
        X = X_;
        Y = Y_;
        direction = direction_;
        vertical = vertical_;
        inner = inner_;
        type = type_;
        xBiasApplied = false;
        yBiasApplied = false;
    }
}

/// <summary>
/// Represents rounding information for geometric elements.
/// </summary>
public class MyRound
{
    /// <summary>Index identifier for this rounding element</summary>
    public int index { get; set; }
    
    /// <summary>Vertical face identifier</summary>
    public int verFace { get; set; }
    
    /// <summary>Horizontal face identifier</summary>
    public int horFace { get; set; }
    
    /// <summary>Maximum radius for this rounding</summary>
    public double MaxRadius { get; set; }
    
    /// <summary>Direction type for the rounding (inner or exterior)</summary>
    public typeRound direction { get; set; }
}