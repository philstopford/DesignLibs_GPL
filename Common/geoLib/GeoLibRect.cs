using System;

namespace geoLib
{
    public class GeoLibRectangle
    {
        public Int32 X { get; set; }
        public Int32 Y { get; set; }
        public Int32 Left { get; set; }
        public Int32 Right { get; set; }
        public Int32 Top { get; set; }
        public Int32 Bottom { get; set; }
        public GeoLibPoint Location { get; set; }
        public Int32 Width { get; set; }
        public Int32 Height { get; set; }

        public GeoLibRectangle()
        {
            pGeoLibRectangle();
        }

        void pGeoLibRectangle()
        {

        }

        public GeoLibRectangle(Int32 x, Int32 y, Int32 width, Int32 height)
        {
            pGeoLibRectangle(x, y, width, height);
        }

        void pGeoLibRectangle(Int32 x, Int32 y, Int32 width, Int32 height)
        {
            X = x;
            Y = y;
            Width = width;
            Height = height;
            setOtherProps();
        }

        public GeoLibRectangle(GeoLibRectangle source)
        {
            pGeoLibRectangle(source);
        }

        void pGeoLibRectangle(GeoLibRectangle source)
        {
            X = source.X;
            Y = source.Y;
            Width = source.Width;
            Height = source.Height;
            setOtherProps();
        }

        void setOtherProps()
        {
            Left = X;
            Top = Y;
            Location = new GeoLibPoint(X, Y);
            Bottom = Top + Height;
            Right = Left + Width;
        }

        public void Offset(GeoLibPoint offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPoint offset)
        {
            X += offset.X;
            Y += offset.Y;
            setOtherProps();
        }
    }

    public class GeoLibRectangleF
    {
        public double X { get; set; }
        public double Y { get; set; }
        public double Left { get; set; }
        public double Right { get; set; }
        public double Top { get; set; }
        public double Bottom { get; set; }
        public GeoLibPointF Location { get; set; }
        public double Width { get; set; }
        public double Height { get; set; }

        public GeoLibRectangleF()
        {
            pGeoLibRectangleF();
        }

        void pGeoLibRectangleF()
        {

        }

        public GeoLibRectangleF(Int32 x, Int32 y, Int32 width, Int32 height)
        {
            pGeoLibRectangleF(x, y, width, height);
        }

        void pGeoLibRectangleF(Int32 x, Int32 y, Int32 width, Int32 height)
        {
            X = x;
            Y = y;
            Width = width;
            Height = height;
            setOtherProps();
        }

        public GeoLibRectangleF(float x, float y, float width, float height)
        {
            pGeoLibRectangleF(x, y, width, height);
        }

        void pGeoLibRectangleF(float x, float y, float width, float height)
        {
            X = x;
            Y = y;
            Width = width;
            Height = height;
            setOtherProps();
        }

        public GeoLibRectangleF(double x, double y, double width, double height)
        {
            pGeoLibRectangleF(x, y, width, height);
        }

        void pGeoLibRectangleF(double x, double y, double width, double height)
        {
            X = (float)x;
            Y = (float)y;
            Width = (float)width;
            Height = (float)height;
            setOtherProps();
        }

        public GeoLibRectangleF(GeoLibRectangle source)
        {
            pGeoLibRectangleF(source);
        }

        void pGeoLibRectangleF(GeoLibRectangle source)
        {
            X = source.X;
            Y = source.Y;
            Width = source.Width;
            Height = source.Height;
            setOtherProps();
        }

        public GeoLibRectangleF(GeoLibRectangleF source)
        {
            pGeoLibRectangleF(source);
        }

        void pGeoLibRectangleF(GeoLibRectangleF source)
        {
            X = source.X;
            Y = source.Y;
            Width = source.Width;
            Height = source.Height;
            setOtherProps();
        }

        void setOtherProps()
        {
            Left = X;
            Top = Y;
            Location = new GeoLibPointF(X, Y);
            Bottom = Top + Height;
            Right = Left + Width;
        }

        public void Offset(GeoLibPoint offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPoint offset)
        {
            X += offset.X;
            Y += offset.Y;
            setOtherProps();
        }

        public void Offset(GeoLibPointF offset)
        {
            pOffset(offset);
        }

        void pOffset(GeoLibPointF offset)
        {
            X += offset.X;
            Y += offset.Y;
            setOtherProps();
        }
    }
}
