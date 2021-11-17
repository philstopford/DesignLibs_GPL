using Eto.Drawing;

namespace etoViewport;

public class ovp_Poly
{
    public PointF[] poly { get; set; }
    public Color color { get; set; }
    public float alpha { get; set; }

    public ovp_Poly(PointF[] geometry, Color geoColor)
    {
        poly = geometry;
        color = geoColor;
        alpha = geoColor.A;
    }

    public ovp_Poly(PointF[] geometry, Color geoColor, float alpha_)
    {
        poly = geometry;
        color = geoColor;
        alpha = alpha_;
    }
}