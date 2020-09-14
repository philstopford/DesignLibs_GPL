using geoLib;

namespace geoCoreLib
{
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

        void pGCStrans()
        {
            reset();
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
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

        void pToggleMirror_x()
        {
            mirror_x = !mirror_x;
            matrix.Scale(1, -1);
        }

        public void setMirror_x()
        {
            pSetMirror_x();
        }

        void pSetMirror_x()
        {
            if (!mirror_x)
            {
                toggleMirror_x();
            }
        }

        public void clearMirror_x()
        {
            pClearMirror_x();
        }

        void pClearMirror_x()
        {
            if (mirror_x)
            {
                toggleMirror_x();
            }
        }

        public void rotate(double a)
        {
            pRotate(a);
        }

        void pRotate(double a)
        {
            if (mirror_x)
            { angle += a; }
            else { angle -= a; }
            if (angle >= 360)
            {
                angle -= 360;
            }
            if (angle < 0)
            {
                angle += 360;
            }
            matrix.Rotate((float)a);
        }

        public void translate(int x, int y)
        {
            pTranslate(x, y);
        }

        void pTranslate(int x, int y)
        {
            matrix.Translate(x, y);
        }

        public void translate(double x, double y)
        {
            pTranslate(x, y);
        }

        void pTranslate(double x, double y)
        {
            matrix.Translate(x, y);
        }

        public void scale(double d)
        {
            pScale(d);
        }

        void pScale(double d)
        {
            mag *= d;
            matrix.Scale((float)d, (float)d);
        }

        public void scale(double x, double y)
        {
            pScale(x, y);
        }

        void pScale(double x, double y)
        {
            mag *= (x + y) / 2.0;
            matrix.Scale((float)x, (float)y);
        }
    }
}
