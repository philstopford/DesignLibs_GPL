using System;

namespace geoLib
{
    public class GeoLibPointF
    {
        public double X { get; set; }
        public double Y { get; set; }

        // Special value for tracking.
        public Int32 tag { get; set; }

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

        public GeoLibPointF(double X_, double Y_)
        {
            pGeoLibPointF(X_, Y_);
        }

        void pGeoLibPointF(double X_, double Y_)
        {
            X = X_;
            Y = Y_;
        }

        public GeoLibPointF(float X_, float Y_)
        {
            pGeoLibPointF(X_, Y_);
        }

        void pGeoLibPointF(float X_, float Y_)
        {
            X = X_;
            Y = Y_;
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

        public Int32 tag { get; set; }

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

        public GeoLibPoint(double X_, double Y_)
        {
            pGeoLibPoint(X_, Y_);
        }

        void pGeoLibPoint(double X_, double Y_)
        {
            X = (Int32)X_;
            Y = (Int32)Y_;
        }

        public GeoLibPoint(float X_, float Y_)
        {
            pGeoLibPoint(X_, Y_);
        }

        void pGeoLibPoint(float X_, float Y_)
        {
            X = (Int32)X_;
            Y = (Int32)Y_;
        }

        public GeoLibPoint(Int32 X_, Int32 Y_)
        {
            pGeoLibPoint(X_, Y_);
        }

        void pGeoLibPoint(Int32 X_, Int32 Y_)
        {
            X = X_;
            Y = Y_;
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
