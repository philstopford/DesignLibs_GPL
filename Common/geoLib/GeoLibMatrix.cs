using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoLib;

public class GeoLibVector3
{
    public double x { get; set; }
    public double y { get; set; }
    public double z { get; set; }

    public GeoLibVector3(GeoLibVector3 source)
    {
        pGeoLibVector3(source);
    }

    private void pGeoLibVector3(GeoLibVector3 source)
    {
        x = source.x;
        y = source.y;
        z = source.z;
    }

    public GeoLibVector3(int X, int Y, int Z)
    {
        pGeoLibVector3(X, Y, Z);
    }

    private void pGeoLibVector3(int X, int Y, int Z)
    {
        x = X;
        y = Y;
        z = Z;
    }

    public GeoLibVector3(double X, double Y, double Z)
    {
        pGeoLibVector3(X, Y, Z);
    }

    private void pGeoLibVector3(double X, double Y, double Z)
    {
        x = X;
        y = Y;
        z = Z;
    }
}

public class GeoLibMatrix
{
    public double[] m { get; set; }

    public GeoLibMatrix(float m11, float m12, float m21, float m22, float dx, float dy)
    {
        pGeoLibMatrix(m11, m12, m21, m22, dx, dy);
    }

    private void pGeoLibMatrix(float m11, float m12, float m21, float m22, float dx, float dy)
    {
        m = new double[6];
        m[0] = m11;
        m[1] = m12;
        m[2] = m21;
        m[3] = m22;
        m[4] = dx;
        m[5] = dy;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Rotate(double ang)
    {
        pRotate(ang);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void pRotate(double ang)
    {
        double m11 = m[0];
        double m22 = m[3];

        // Cache sin/cos calculations and convert to radians once
        double angleRad = ang * Math.PI / 180.0;
        double cosAngle = Math.Cos(angleRad);
        double sinAngle = Math.Sin(angleRad);

        m[0] = cosAngle * m11;
        m[1] = sinAngle * m22;
        m[2] = -sinAngle * m11;
        m[3] = cosAngle * m22;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public void Scale(float sx, float sy)
    {
        pScale(sx, sy);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void pScale(float sx, float sy)
    {
        m[0] *= sx;
        m[3] *= sy;
    }

    public void Translate(int x, int y)
    {
        pTranslate(x, y);
    }

    private void pTranslate(int x, int y)
    {
        m[4] += x;
        m[5] += y;
    }

    public void Translate(double x, double y)
    {
        pTranslate(x, y);
    }

    private void pTranslate(double x, double y)
    {
        m[4] += x;
        m[5] += y;
    }

    public void TransformPoints(Path64 source)
    {
        pTransformPoints(source);
    }

    private void pTransformPoints(Path64 source)
    {
        double m11 = m[0];
        double m12 = m[1];
        double m21 = m[2];
        double m22 = m[3];
        double dx = m[4];
        double dy = m[5];

        int sLength = source.Count;
#if !GEOLIBSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                double x1 = m11 * source[pt].X + m21 * source[pt].Y + dx;
                double y1 = m12 * source[pt].X + m22 * source[pt].Y + dy;
                source[pt] = new Point64(x1, y1);
            }
#if !GEOLIBSINGLETHREADED
        );
#endif
    }

    public void TransformPoints(PathD source)
    {
        pTransformPoints(source);
    }

    private void pTransformPoints(PathD source)
    {
        double m11 = m[0];
        double m12 = m[1];
        double m21 = m[2];
        double m22 = m[3];
        double dx = m[4];
        double dy = m[5];

        int sLength = source.Count;
#if !GEOLIBSINGLETHREADED
        Parallel.For(0, sLength, pt => 
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                double x1 = m11 * source[pt].x + m21 * source[pt].y + dx;
                double y1 = m12 * source[pt].x + m22 * source[pt].y + dy;
                source[pt] = new PointD(x1, y1);
            }
#if !GEOLIBSINGLETHREADED
        );
#endif
    }
}