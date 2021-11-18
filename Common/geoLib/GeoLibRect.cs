namespace geoLib;

public class GeoLibRectangle
{
    public int X { get; set; }
    public int Y { get; set; }
    public int Left { get; set; }
    public int Right { get; set; }
    public int Top { get; set; }
    public int Bottom { get; set; }
    public GeoLibPoint Location { get; set; }
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
        X = x;
        Y = y;
        Width = width;
        Height = height;
        setOtherProps();
    }

    public GeoLibRectangle(GeoLibRectangle source)
    {
        pGeoLibRectangle(source);
    }

    private void pGeoLibRectangle(GeoLibRectangle source)
    {
        X = source.X;
        Y = source.Y;
        Width = source.Width;
        Height = source.Height;
        setOtherProps();
    }

    private void setOtherProps()
    {
        Left = X;
        Top = Y;
        Location = new GeoLibPoint(X, Y);
        Bottom = Top + Height;
        Right = Left + Width;
    }

    public void Offset(GeoLibPoint offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint offset)
    {
        X += offset.X;
        Y += offset.Y;
        setOtherProps();
    }
}

public class GeoLibRectangleF
{
    public double X { get; set; }
    public double Y { get; set; }
    public double Left { get; set; }
    public double Right { get; set; }
    public double Top { get; set; }
    public double Bottom { get; set; }
    public GeoLibPointF Location { get; set; }
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
        X = x;
        Y = y;
        Width = width;
        Height = height;
        setOtherProps();
    }

    public GeoLibRectangleF(float x, float y, float width, float height)
    {
        pGeoLibRectangleF(x, y, width, height);
    }

    private void pGeoLibRectangleF(float x, float y, float width, float height)
    {
        X = x;
        Y = y;
        Width = width;
        Height = height;
        setOtherProps();
    }

    public GeoLibRectangleF(double x, double y, double width, double height)
    {
        pGeoLibRectangleF(x, y, width, height);
    }

    private void pGeoLibRectangleF(double x, double y, double width, double height)
    {
        X = (float)x;
        Y = (float)y;
        Width = (float)width;
        Height = (float)height;
        setOtherProps();
    }

    public GeoLibRectangleF(GeoLibRectangle source)
    {
        pGeoLibRectangleF(source);
    }

    private void pGeoLibRectangleF(GeoLibRectangle source)
    {
        X = source.X;
        Y = source.Y;
        Width = source.Width;
        Height = source.Height;
        setOtherProps();
    }

    public GeoLibRectangleF(GeoLibRectangleF source)
    {
        pGeoLibRectangleF(source);
    }

    private void pGeoLibRectangleF(GeoLibRectangleF source)
    {
        X = source.X;
        Y = source.Y;
        Width = source.Width;
        Height = source.Height;
        setOtherProps();
    }

    private void setOtherProps()
    {
        Left = X;
        Top = Y;
        Location = new GeoLibPointF(X, Y);
        Bottom = Top + Height;
        Right = Left + Width;
    }

    public void Offset(GeoLibPoint offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint offset)
    {
        X += offset.X;
        Y += offset.Y;
        setOtherProps();
    }

    public void Offset(GeoLibPointF offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPointF offset)
    {
        X += offset.X;
        Y += offset.Y;
        setOtherProps();
    }
}