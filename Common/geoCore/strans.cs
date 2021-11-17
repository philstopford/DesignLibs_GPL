using geoLib;

namespace geoCoreLib;

public class GCStrans
{
    public GeoLibMatrix matrix { get; set; }
    public double mag { get; set; }
    public double angle { get; set; }
    public bool mirror_x { get; set; }

    public GCStrans()
    {
        pGCStrans();
    }

    private void pGCStrans()
    {
        reset();
    }

    public void reset()
    {
        pReset();
    }

    private void pReset()
    {
        mag = 1;
        angle = 0;
        mirror_x = false;
        matrix = new GeoLibMatrix(1, 0, 0, 1, 0, 0);
    }

    public void toggleMirror_x()
    {
        pToggleMirror_x();
    }

    private void pToggleMirror_x()
    {
        mirror_x = !mirror_x;
        matrix.Scale(1, -1);
    }

    public void setMirror_x()
    {
        pSetMirror_x();
    }

    private void pSetMirror_x()
    {
        switch (mirror_x)
        {
            case false:
                toggleMirror_x();
                break;
        }
    }

    public void clearMirror_x()
    {
        pClearMirror_x();
    }

    private void pClearMirror_x()
    {
        switch (mirror_x)
        {
            case true:
                toggleMirror_x();
                break;
        }
    }

    public void rotate(double a)
    {
        pRotate(a);
    }

    private void pRotate(double a)
    {
        switch (mirror_x)
        {
            case true:
                angle += a;
                break;
            default:
                angle -= a;
                break;
        }
        switch (angle)
        {
            case >= 360:
                angle -= 360;
                break;
        }
        switch (angle)
        {
            case < 0:
                angle += 360;
                break;
        }
        matrix.Rotate((float)a);
    }

    public void translate(int x, int y)
    {
        pTranslate(x, y);
    }

    private void pTranslate(int x, int y)
    {
        matrix.Translate(x, y);
    }

    public void translate(double x, double y)
    {
        pTranslate(x, y);
    }

    private void pTranslate(double x, double y)
    {
        matrix.Translate(x, y);
    }

    public void scale(double d)
    {
        pScale(d);
    }

    private void pScale(double d)
    {
        mag *= d;
        matrix.Scale((float)d, (float)d);
    }

    public void scale(double x, double y)
    {
        pScale(x, y);
    }

    private void pScale(double x, double y)
    {
        mag *= (x + y) / 2.0;
        matrix.Scale((float)x, (float)y);
    }
}