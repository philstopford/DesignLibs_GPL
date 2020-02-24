using System;
using System.Threading.Tasks;

namespace geoLib
{
    public class GeoLibVector3
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }

        public GeoLibVector3(GeoLibVector3 source)
        {
            pGeoLibVector3(source);
        }

        void pGeoLibVector3(GeoLibVector3 source)
        {
            x = source.x;
            y = source.y;
            z = source.z;
        }

        public GeoLibVector3(Int32 X, Int32 Y, Int32 Z)
        {
            pGeoLibVector3(X, Y, Z);
        }

        void pGeoLibVector3(Int32 X, Int32 Y, Int32 Z)
        {
            x = X;
            y = Y;
            z = Z;
        }

        public GeoLibVector3(double X, double Y, double Z)
        {
            pGeoLibVector3(X, Y, Z);
        }

        void pGeoLibVector3(double X, double Y, double Z)
        {
            x = X;
            y = Y;
            z = Z;
        }
    }

    public class GeoLibMatrix
    {
        public double[] m { get; set; }

        public GeoLibMatrix()
        {
            pGeoLibMatrix();
        }

        void pGeoLibMatrix()
        {
        }

        public GeoLibMatrix(float m11, float m12, float m21, float m22, float dx, float dy)
        {
            pGeoLibMatrix(m11, m12, m21, m22, dx, dy);
        }

        void pGeoLibMatrix(float m11, float m12, float m21, float m22, float dx, float dy)
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

        void pRotate(double ang)
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

        void pScale(float sx, float sy)
        {
            m[0] *= sx;
            m[3] *= sy;
        }

        public void Translate(Int32 x, Int32 y)
        {
            pTranslate(x, y);
        }

        void pTranslate(Int32 x, Int32 y)
        {
            m[4] += x;
            m[5] += y;
        }

        public void Translate(double x, double y)
        {
            pTranslate(x, y);
        }

        void pTranslate(double x, double y)
        {
            m[4] += x;
            m[5] += y;
        }

        public void TransformPoints(GeoLibPoint[] source)
        {
            pTransformPoints(source);
        }

        void pTransformPoints(GeoLibPoint[] source)
        {
            var m11 = m[0];
            var m12 = m[1];
            var m21 = m[2];
            var m22 = m[3];
            var dx = m[4];
            var dy = m[5];

            int sLength = source.Length;
            Parallel.For(0, sLength, (pt) => // (int pt = 0; pt < source.Length; pt++)
            {
                double x1 = m11 * source[pt].X + m21 * source[pt].Y + dx;
                double y1 = m12 * source[pt].X + m22 * source[pt].Y + dy;
                source[pt] = new GeoLibPoint(x1, y1);
            });
        }

        public void TransformPoints(GeoLibPointF[] source)
        {
            pTransformPoints(source);
        }

        void pTransformPoints(GeoLibPointF[] source)
        {
            var m11 = m[0];
            var m12 = m[1];
            var m21 = m[2];
            var m22 = m[3];
            var dx = m[4];
            var dy = m[5];

            int sLength = source.Length;
            Parallel.For(0, sLength, (pt) => // (int pt = 0; pt < source.Length; pt++)
            {
                double x1 = m11 * source[pt].X + m21 * source[pt].Y + dx;
                double y1 = m12 * source[pt].X + m22 * source[pt].Y + dy;
                source[pt] = new GeoLibPointF(x1, y1);
            });
        }
    }
}
