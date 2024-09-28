using Clipper2Lib;

namespace geoLib;

public class GeoLibRectangle
{
    public Point64 Location { get; set; }
    public int Width { get; set; }
    public int Height { get; set; }

    public GeoLibRectangle()
    {
        pGeoLibRectangle();
    }

    private void pGeoLibRectangle()
    {

    }

    public GeoLibRectangle(int x, int y, int width, int height)
    {
        pGeoLibRectangle(x, y, width, height);
    }

    private void pGeoLibRectangle(int x, int y, int width, int height)
    {
        Location = new Point64(x, y);
        Width = width;
        Height = height;
    }

    public GeoLibRectangle(GeoLibRectangle source)
    {
        pGeoLibRectangle(source);
    }

    private void pGeoLibRectangle(GeoLibRectangle source)
    {
        Location = new Point64(source.Location);
        Width = source.Width;
        Height = source.Height;
    }
    
    public void Offset(Point64 offset)
    {
        pOffset(offset);
    }

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