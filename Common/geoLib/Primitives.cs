using System;

namespace geoLib
{
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

        void pMyVertex(double X, double Y, typeDirection direction, bool vertical, bool inner, typeVertex type)
        {
            this.X = X;
            this.Y = Y;
            this.direction = direction;
            this.vertical = vertical;
            this.inner = inner;
            this.type = type;
            xBiasApplied = false;
            yBiasApplied = false;
        }
    }

    public class MyRound
    {
        public Int32 index { get; set; }
        public Int32 verFace { get; set; }
        public Int32 horFace { get; set; }
        public double MaxRadius { get; set; }
        public typeRound direction { get; set; }

        public MyRound()
        {
            pMyRound();
        }

        void pMyRound()
        {
            // nothing to do here.
        }
    }
}
