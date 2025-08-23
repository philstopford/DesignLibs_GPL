using Clipper2Lib;

namespace geoLib;

/// <summary>
/// Rectangle implementation to avoid platform-specific quirks of System.Drawing.
/// Provides basic rectangle operations for multi-platform compatibility.
/// </summary>
public class GeoLibRectangle
{
    /// <summary>Location of the rectangle's top-left corner</summary>
    public Point64 Location { get; set; }
    
    /// <summary>Width of the rectangle</summary>
    public int Width { get; set; }
    
    /// <summary>Height of the rectangle</summary>
    public int Height { get; set; }

    /// <summary>
    /// Creates a new rectangle at origin with zero dimensions.
    /// </summary>
    public GeoLibRectangle()
    {
        pGeoLibRectangle();
    }

    /// <summary>Internal initialization for default constructor</summary>
    private void pGeoLibRectangle()
    {

    }

    /// <summary>
    /// Creates a new rectangle with specified location and dimensions.
    /// </summary>
    /// <param name="x">X coordinate of top-left corner</param>
    /// <param name="y">Y coordinate of top-left corner</param>
    /// <param name="width">Width of rectangle</param>
    /// <param name="height">Height of rectangle</param>
    public GeoLibRectangle(int x, int y, int width, int height)
    {
        pGeoLibRectangle(x, y, width, height);
    }

    /// <summary>Internal initialization with parameters</summary>
    private void pGeoLibRectangle(int x, int y, int width, int height)
    {
        Location = new Point64(x, y);
        Width = width;
        Height = height;
    }

    /// <summary>
    /// Creates a copy of an existing rectangle.
    /// </summary>
    /// <param name="source">Rectangle to copy from</param>
    public GeoLibRectangle(GeoLibRectangle source)
    {
        pGeoLibRectangle(source);
    }

    /// <summary>Internal initialization from source rectangle</summary>
    private void pGeoLibRectangle(GeoLibRectangle source)
    {
        Location = new Point64(source.Location);
        Width = source.Width;
        Height = source.Height;
    }
    
    /// <summary>
    /// Moves the rectangle by the specified offset.
    /// </summary>
    /// <param name="offset">Amount to move the rectangle</param>
    public void Offset(Point64 offset)
    {
        pOffset(offset);
    }

    /// <summary>Internal implementation of offset operation</summary>
    private void pOffset(Point64 offset)
    {
        Location = new Point64(Location.X + offset.X, Location.Y + offset.Y);
    }
}

public class GeoLibRectangleF
{
    public PointD Location { get; set; }
    public double Width { get; set; }
    public double Height { get; set; }

    public GeoLibRectangleF()
    {
        pGeoLibRectangleF();
    }

    private void pGeoLibRectangleF()
    {

    }

    public GeoLibRectangleF(int x, int y, int width, int height)
    {
        pGeoLibRectangleF(x, y, width, height);
    }

    private void pGeoLibRectangleF(int x, int y, int width, int height)
    {
        Location = new PointD(x, y);
        Width = width;
        Height = height;
    }

    public GeoLibRectangleF(float x, float y, float width, float height)
    {
        pGeoLibRectangleF(x, y, width, height);
    }

    private void pGeoLibRectangleF(float x, float y, float width, float height)
    {
        Location = new PointD(x, y);
        Width = width;
        Height = height;
    }

    public GeoLibRectangleF(double x, double y, double width, double height)
    {
        pGeoLibRectangleF(x, y, width, height);
    }

    private void pGeoLibRectangleF(double x, double y, double width, double height)
    {
        Location = new PointD(x, y);
        Width = (float)width;
        Height = (float)height;
    }

    public GeoLibRectangleF(GeoLibRectangle source)
    {
        pGeoLibRectangleF(source);
    }

    private void pGeoLibRectangleF(GeoLibRectangle source)
    {
        Location = new PointD(source.Location);
        Width = source.Width;
        Height = source.Height;
    }

    public GeoLibRectangleF(GeoLibRectangleF source)
    {
        pGeoLibRectangleF(source);
    }

    private void pGeoLibRectangleF(GeoLibRectangleF source)
    {
        Location = new PointD(source.Location);
        Width = source.Width;
        Height = source.Height;
    }
    
    public void Offset(Point64 offset)
    {
        pOffset(offset);
    }

    private void pOffset(Point64 offset)
    {
        Location = new PointD(Location.x + offset.X, Location.y + offset.Y);
    }

    public void Offset(PointD offset)
    {
        pOffset(offset);
    }

    private void pOffset(PointD offset)
    {
        Location = new PointD(Location.x + offset.x, Location.y + offset.y);
    }
}