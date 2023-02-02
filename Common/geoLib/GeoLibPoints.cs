using Clipper2Lib;

namespace geoLib;

public class GeoLibPointF_
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

    public GeoLibPointF_()
    {
        pGeoLibPointF();
    }

    private void pGeoLibPointF()
    {
        internal_pt = new();
    }

    public GeoLibPointF_(GeoLibPointF_ sourcePoint)
    {
        pGeoLibPointF(sourcePoint);
    }

    private void pGeoLibPointF(GeoLibPointF_ sourcePoint)
    {
        internal_pt = new();
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPointF_(GeoLibPoint_ sourcePoint)
    {
        pGeoLibPointF(sourcePoint);
    }

    private void pGeoLibPointF(GeoLibPoint_ sourcePoint)
    {
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPointF_(double X_, double Y_)
    {
        pGeoLibPointF(X_, Y_);
    }

    private void pGeoLibPointF(double X_, double Y_)
    {
        internal_pt = new();
        X = X_;
        Y = Y_;
    }

    public GeoLibPointF_(float X_, float Y_)
    {
        pGeoLibPointF((double)X_, (double)Y_);
    }
    
    public void Offset(GeoLibPoint_ offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint_ offset)
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

    public void Offset(GeoLibPointF_ offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPointF_ offset)
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

public class GeoLibPoint_
{
    private Point64 internal_pt;

    public long X
    {
        get => internal_pt.X;

        set => internal_pt.X = value;
    }
    
    public long Y
    {
        get => internal_pt.Y;

        set => internal_pt.Y = value;
    }
    
    public int tag { get; set; }

    public GeoLibPoint_()
    {
        pGeoLibPoint();
    }

    private void pGeoLibPoint()
    {

    }

    public GeoLibPoint_(GeoLibPoint_ sourcePoint)
    {
        pGeoLibPoint(sourcePoint);
    }

    private void pGeoLibPoint(GeoLibPoint_ sourcePoint)
    {
        X = sourcePoint.X;
        Y = sourcePoint.Y;
    }

    public GeoLibPoint_(double X_, double Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(double X_, double Y_)
    {
        X = (int)X_;
        Y = (int)Y_;
    }

    public GeoLibPoint_(float X_, float Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(float X_, float Y_)
    {
        X = (int)X_;
        Y = (int)Y_;
    }

    public GeoLibPoint_(int X_, int Y_)
    {
        pGeoLibPoint(X_, Y_);
    }

    private void pGeoLibPoint(int X_, int Y_)
    {
        X = X_;
        Y = Y_;
    }

    public void Offset(GeoLibPoint_ offset)
    {
        pOffset(offset);
    }

    private void pOffset(GeoLibPoint_ offset)
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