namespace geoLib;

public enum typeDirection { left1, right1, up1, down1, tilt1 };
public enum typeVertex { corner, center };
public enum typeRound { inner, exter };

public class MyVertex : GeoLibPointF
{
    public typeDirection direction { get; set; }
    public bool vertical { get; set; } // denotes a vertex attached to a vertical edge
    public bool inner { get; set; } // denotes an inner vertex
    public typeVertex type { get; set; } // this defines whether the point is a vertex or a center point on an edge. 

    // Track whether we applied bias to this vertex already.
    public bool xBiasApplied { get; set; }
    public bool yBiasApplied { get; set; }

    public MyVertex(double X, double Y, typeDirection direction, bool vertical, bool inner, typeVertex type)
    {
        pMyVertex(X, Y, direction, vertical, inner, type);
    }

    public MyVertex(MyVertex source)
    {
        pMyVertex(source.X, source.Y, source.direction, source.vertical, source.inner, source.type);
    }

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

public class MyRound
{
    public int index { get; set; }
    public int verFace { get; set; }
    public int horFace { get; set; }
    public double MaxRadius { get; set; }
    public typeRound direction { get; set; }

    public MyRound()
    {
        pMyRound();
    }

    private void pMyRound()
    {
        // nothing to do here.
    }
}