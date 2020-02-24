using System;

namespace geoLib
{
    public class GeoLibPointF
    {
        public double X { get; set; }
        public double Y { get; set; }

        public GeoLibPointF()
        {
            pGeoLibPointF();
        }

        void pGeoLibPointF()
        {

        }

        public GeoLibPointF(GeoLibPointF sourcePoint)
        {
            pGeoLibPointF(sourcePoint);
        }

        void pGeoLibPointF(GeoLibPointF sourcePoint)
        {
            X = sourcePoint.X;
            Y = sourcePoint.Y;
        }

        public GeoLibPointF(GeoLibPoint sourcePoint)
        {
            pGeoLibPointF(sourcePoint);
        }

        void pGeoLibPointF(GeoLibPoint sourcePoint)
        {
            X = sourcePoint.X;
            Y = sourcePoint.Y;
        }

        public GeoLibPointF(double X, double Y)
        {
            pGeoLibPointF(X, Y);
        }

        void pGeoLibPointF(double X, double Y)
        {
            this.X = X;
            this.Y = Y;
        }

        public GeoLibPointF(float X, float Y)
        {
            pGeoLibPointF(X, Y);
        }

        void pGeoLibPointF(float X, float Y)
        {
            this.X = X;
            this.Y = Y;
        }

        public void Offset(GeoLibPoint offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPoint offset)
        {
            X += offset.X;
            Y += offset.Y;
        }

        public void Offset(Int32 x, Int32 y)
        {
            pOffset(x, y);
        }

        void pOffset(Int32 x, Int32 y)
        {
            X += x;
            Y += y;
        }

        public void Offset(GeoLibPointF offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPointF offset)
        {
            X += offset.X;
            Y += offset.Y;
        }

        public void Offset(double x, double y)
        {
            pOffset(x, y);
        }

        void pOffset(double x, double y)
        {
            X += x;
            Y += y;
        }
    }

    public class GeoLibPoint
    {
        public Int32 X { get; set; }
        public Int32 Y { get; set; }

        public GeoLibPoint()
        {
            pGeoLibPoint();
        }

        void pGeoLibPoint()
        {

        }

        public GeoLibPoint(GeoLibPoint sourcePoint)
        {
            pGeoLibPoint(sourcePoint);
        }

        void pGeoLibPoint(GeoLibPoint sourcePoint)
        {
            X = sourcePoint.X;
            Y = sourcePoint.Y;
        }

        public GeoLibPoint(double X, double Y)
        {
            pGeoLibPoint(X, Y);
        }

        void pGeoLibPoint(double X, double Y)
        {
            this.X = (Int32)X;
            this.Y = (Int32)Y;
        }

        public GeoLibPoint(float X, float Y)
        {
            pGeoLibPoint(X, Y);
        }

        void pGeoLibPoint(float X, float Y)
        {
            this.X = (Int32)X;
            this.Y = (Int32)Y;
        }

        public GeoLibPoint(Int32 X, Int32 Y)
        {
            pGeoLibPoint(X, Y);
        }

        void pGeoLibPoint(Int32 X, Int32 Y)
        {
            this.X = X;
            this.Y = Y;
        }

        public void Offset(GeoLibPoint offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPoint offset)
        {
            X += offset.X;
            Y += offset.Y;
        }

        public void Offset(Int32 x, Int32 y)
        {
            pOffset(x, y);
        }

        void pOffset(Int32 x, Int32 y)
        {
            X += x;
            Y += y;
        }
    }
}
