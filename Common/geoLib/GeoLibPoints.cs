using Clipper2Lib;

namespace geoLib;

public class GeoLibPointF
{
    private PointD internal_pt;

    public double X
    {
        get => internal_pt.x;

        set => internal_pt.x = value;
    }
    
    public double Y
    {
        get => internal_pt.y;

        set => internal_pt.y = value;
    }

    // Special value for tracking.
    public int tag { get; set; }

    public GeoLibPointF()
    {
        pGeoLibPointF();
    }

    private void pGeoLibPointF()
    {
        internal_pt = new();
    }

    public GeoLibPointF(GeoLibPointF sourcePoint)
    {
        pGeoLibPointF(sourcePoint);
    }

    private void pGeoLibPointF(GeoLibPointF sourcePoint)
    {
        internal_pt = new();
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPointF(GeoLibPoint sourcePoint)
    {
        pGeoLibPointF(sourcePoint);
    }

    private void pGeoLibPointF(GeoLibPoint sourcePoint)
    {
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPointF(double X_, double Y_)
    {
        pGeoLibPointF(X_, Y_);
    }

    private void pGeoLibPointF(double X_, double Y_)
    {
        internal_pt = new();
        X = X_;
        Y = Y_;
    }

    public GeoLibPointF(float X_, float Y_)
    {
        pGeoLibPointF((double)X_, (double)Y_);
    }
    
    public void Offset(GeoLibPoint offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint offset)
    {
        X += offset.X;
        Y += offset.Y;
    }

    public void Offset(int x, int y)
    {
        pOffset(x, y);
    }

    private void pOffset(int x, int y)
    {
        X += x;
        Y += y;
    }

    public void Offset(GeoLibPointF offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPointF offset)
    {
        X += offset.X;
        Y += offset.Y;
    }

    public void Offset(double x, double y)
    {
        pOffset(x, y);
    }

    private void pOffset(double x, double y)
    {
        X += x;
        Y += y;
    }
}

public class GeoLibPoint
{
    public int X { get; set; }
    public int Y { get; set; }

    public int tag { get; set; }

    public GeoLibPoint()
    {
        pGeoLibPoint();
    }

    private void pGeoLibPoint()
    {

    }

    public GeoLibPoint(GeoLibPoint sourcePoint)
    {
        pGeoLibPoint(sourcePoint);
    }

    private void pGeoLibPoint(GeoLibPoint sourcePoint)
    {
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPoint(double X_, double Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(double X_, double Y_)
    {
        X = (int)X_;
        Y = (int)Y_;
    }

    public GeoLibPoint(float X_, float Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(float X_, float Y_)
    {
        X = (int)X_;
        Y = (int)Y_;
    }

    public GeoLibPoint(int X_, int Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(int X_, int Y_)
    {
        X = X_;
        Y = Y_;
    }

    public void Offset(GeoLibPoint offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint offset)
    {
        X += offset.X;
        Y += offset.Y;
    }

    public void Offset(int x, int y)
    {
        pOffset(x, y);
    }

    private void pOffset(int x, int y)
    {
        X += x;
        Y += y;
    }
}