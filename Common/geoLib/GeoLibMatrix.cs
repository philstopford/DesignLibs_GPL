using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoLib;

/// <summary>
/// 3D vector implementation for geometric calculations.
/// </summary>
public class GeoLibVector3
{
    /// <summary>X component of the vector</summary>
    public double x { get; set; }
    
    /// <summary>Y component of the vector</summary>
    public double y { get; set; }
    
    /// <summary>Z component of the vector</summary>
    public double z { get; set; }

    /// <summary>
    /// Creates a copy of an existing vector.
    /// </summary>
    /// <param name="source">Vector to copy from</param>
    public GeoLibVector3(GeoLibVector3 source)
    {
        pGeoLibVector3(source);
    }

    /// <summary>Internal copy initialization</summary>
    private void pGeoLibVector3(GeoLibVector3 source)
    {
        x = source.x;
        y = source.y;
        z = source.z;
    }

    /// <summary>
    /// Creates a new vector with integer coordinates.
    /// </summary>
    /// <param name="X">X component</param>
    /// <param name="Y">Y component</param>
    /// <param name="Z">Z component</param>
    public GeoLibVector3(int X, int Y, int Z)
    {
        pGeoLibVector3(X, Y, Z);
    }

    /// <summary>Internal initialization with integer values</summary>
    private void pGeoLibVector3(int X, int Y, int Z)
    {
        x = X;
        y = Y;
        z = Z;
    }

    /// <summary>
    /// Creates a new vector with floating-point coordinates.
    /// </summary>
    /// <param name="X">X component</param>
    /// <param name="Y">Y component</param>
    /// <param name="Z">Z component</param>
    public GeoLibVector3(double X, double Y, double Z)
    {
        pGeoLibVector3(X, Y, Z);
    }

    /// <summary>Internal initialization with double values</summary>
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

    public void Rotate(double ang)
    {
        pRotate(ang);
    }

    private void pRotate(double ang)
    {
        double m11 = m[0];
        double m22 = m[3];

        m[0] = Math.Cos(ang * Math.PI / 180) * m11;
        m[1] = Math.Sin(ang * Math.PI / 180) * m22;
        m[2] = -Math.Sin(ang * Math.PI / 180) * m11;
        m[3] = Math.Cos(ang * Math.PI / 180) * m22;
    }

    public void Scale(float sx, float sy)
    {
        pScale(sx, sy);
    }

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