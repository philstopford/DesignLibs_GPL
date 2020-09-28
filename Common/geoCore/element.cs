using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;

namespace geoCoreLib
{
    public class GCElement
    {
        public Int32 layer_nr { get; set; }
        public Int32 datatype_nr { get; set; }
        public bool select { get; set; }

        public double angle(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            return pAngle(p1, p2, p3);
        }

        double pAngle(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            double a1 = 0, a2 = 0;
            GeoLibPointF dif1, dif2;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            dif2 = new GeoLibPointF(p3.X - p2.X, p3.Y - p2.Y);
            if (dif1.X != 0)
            {
                a1 = Math.Atan(dif1.Y / dif1.X) / 2 / Math.PI * 360;
                if (dif1.X < 0) { a1 -= 180; }
                if (a1 < -180) { a1 += 360; }
            }
            else { if (dif1.Y > 0) { a1 = 90; } else { a1 = -90; } }
            if (dif2.X != 0)
            {
                a2 = Math.Atan(dif2.Y / dif2.X) / 2 / Math.PI * 360;
                if (dif2.X < 0) { a2 -= 180; }
                if (a2 < -180) { a2 += 360; }
            }
            else { if (dif2.Y > 0) { a2 = 90; } else { a2 = -90; } }
            a1 = a2 - a1;
            if (a1 <= -180) a1 += 360;
            if (a1 >= 180) a1 -= 360;
            return (a1);
        }

        public double angle(GeoLibPoint p1, GeoLibPoint p2)
        {
            return pAngle(p1, p2);
        }

        double pAngle(GeoLibPoint p1, GeoLibPoint p2)
        {
            double a1 = 0;
            GeoLibPointF dif1;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            if (dif1.X != 0)
            {
                a1 = Math.Atan(dif1.Y / dif1.X) / 2 / Math.PI * 360;
                if (dif1.X < 0) { a1 -= 180; }
                if (a1 < -180) { a1 += 360; }
            }
            else { if (dif1.Y > 0) { a1 = 90; } else { a1 = -90; } }
            if (a1 <= -180) a1 += 360;
            if (a1 > 180) a1 -= 360;
            return (a1);
        }

        public bool cutPoint2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4, GeoLibPoint pc)
        {
            return pCutPoint2(p1, p2, p3, p4, pc);
        }

        bool pCutPoint2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4, GeoLibPoint pc)
        {
            //between  p1-p2+p3-p4
            double m1, m2, b1, b2;
            GeoLibPointF dif1, dif2;
            int help;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            dif2 = new GeoLibPointF(p4.X - p3.X, p4.Y - p3.Y);
            pc.X = (1 << 30);
            pc.Y = (1 << 30);
            if ((dif1.X != 0) && (dif2.X != 0))
            {
                m1 = dif1.Y / dif1.X;
                m2 = dif2.Y / dif2.X;
                b1 = -m1 * p1.X + p1.Y;
                b2 = -m2 * p3.X + p3.Y;
                if (m1 == m2) { return false; }
                pc.X = (round((b2 - b1) / (m1 - m2)));
                if (dif1.Y == 0) { pc.Y = (p1.Y); }
                else if (dif2.Y == 0) { pc.Y = (p3.Y); }
                else pc.Y = (round(m1 * pc.X + b1));
            }
            if ((dif1.X == 0) && (dif2.X != 0))
            {
                m2 = dif2.Y / dif2.X;
                b2 = -m2 * p3.X + p3.Y;
                pc.X = (p1.X);
                if (dif2.Y != 0) { pc.Y = (round(m2 * pc.X + b2)); }
                else pc.Y = (p3.Y);
            }
            if ((dif1.X != 0) && (dif2.X == 0))
            {
                m1 = dif1.Y / dif1.X;
                b1 = -m1 * p1.X + p1.Y;
                pc.X = (p3.X);
                if (dif2.Y != 0) { pc.Y = (round(m1 * pc.X + b1)); }
                else pc.Y = (p1.Y);
            }
            if ((dif1.X == 0) && (dif2.X == 0)) { return false; }
            if (pc == p1) { return false; }
            if (pc == p2) { return false; }
            if (pc == p3) { return false; }
            if (pc == p4) { return false; }
            if (p1.X > p2.X) { help = p1.X; p1.X = (p2.X); p2.X = (help); }
            if (p3.X > p4.X) { help = p3.X; p3.X = (p4.X); p4.X = (help); }
            if (pc.X < p1.X) { return false; }
            if (pc.X > p2.X) { return false; }
            if (pc.X < p3.X) { return false; }
            if (pc.X > p4.X) { return false; }
            if (p1.Y > p2.Y) { help = p1.Y; p1.Y = (p2.Y); p2.Y = (help); }
            if (p3.Y > p4.Y) { help = p3.Y; p3.Y = (p4.Y); p4.Y = (help); }
            if (pc.Y < p1.Y) { return false; }
            if (pc.Y > p2.Y) { return false; }
            if (pc.Y < p3.Y) { return false; }
            if (pc.Y > p4.Y) { return false; }
            return true;
        }

        public double distance(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            return pDistance(p1, p2, p3);
        }

        double pDistance(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            // >0 if left off line
            double m1, b1, wert;
            GeoLibPoint dif1;
            dif1 = new GeoLibPoint(p2.X - p1.X, p2.Y - p1.Y);
            if ((dif1.X != 0))
            {
                m1 = (double)(dif1.Y) / (dif1.X);
                b1 = -m1 * p1.X + p1.Y;
                wert = m1 * p3.X + b1 - p3.Y;
                if (dif1.Y == 0)
                {
                    if (dif1.X > 0)
                    {
                        return -wert;
                    }
                    else
                    {
                        return wert;
                    }
                }
                if (dif1.X > 0)
                {
                    return -wert;
                }
                if (dif1.X < 0)
                {
                    return +wert;
                }
                else
                {
                    return wert;
                }
            }
            else
            {
                if (dif1.Y > 0)
                {
                    return (p1.X - p3.X);
                }
                else
                {
                    return (p3.X - p1.X);
                }
            }
        }

        public double distance(GeoLibPoint p1, GeoLibPoint p2)
        {
            return pDistance(p1, p2);
        }

        double pDistance(GeoLibPoint p1, GeoLibPoint p2)
        {
            Int32 dx = p1.X - p2.X;
            Int32 dy = p1.Y - p2.Y;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        public bool identical(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            return pIdentical(p1, p2, p3, p4);
        }

        bool pIdentical(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            double m1, m2, b1, b2;
            GeoLibPointF dif1, dif2;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            dif2 = new GeoLibPointF(p4.X - p3.X, p4.Y - p3.Y);
            if ((dif1.Y != 0) && (dif2.Y != 0))
            {
                m1 = dif1.X / dif1.Y;
                m2 = dif2.X / dif2.Y;
                //printf("m %f %f\n",m1,m2);
                if (m1 != m2) { return false; }
                b1 = -m1 * p1.Y + p1.X;
                b2 = -m2 * p3.Y + p3.X;
                //printf("b %f %f\n",b1,b2);
                if (b1 != b2) { return false; }
                return true;
            }
            else if ((dif1.X != 0) && (dif2.X != 0))
            {
                m1 = dif1.Y / dif1.X;
                m2 = dif2.Y / dif2.X;
                //printf("m %f %f\n",m1,m2);
                if (m1 != m2) { return false; }
                b1 = -m1 * p1.X + p1.Y;
                b2 = -m2 * p3.X + p3.Y;
                //printf("b %f %f\n",b1,b2);
                if (b1 != b2) { return false; }
                return true;
            }
            return false;
        }

        public double length(GeoLibPoint p)
        {
            return pLength(p);
        }

        double pLength(GeoLibPoint p)
        {
            return Math.Sqrt(p.X * p.X + p.Y * p.Y);
        }

        public bool parallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            return pParallel(p1, p2, p3, p4);
        }

        bool pParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            double m1, m2;
            GeoLibPoint dif1, dif2;
            dif1 = new GeoLibPoint(p2.X - p1.X, p2.Y - p1.Y);
            dif2 = new GeoLibPoint(p4.X - p3.X, p4.Y - p3.Y);
            if ((dif1.Y != 0) && (dif2.Y != 0))
            {
                m1 = (double)(dif1.X) / (dif1.Y);
                m2 = (double)(dif2.X) / (dif2.Y);
                if (m1 != m2)
                {
                    return false;
                }
                return true;
            }
            else if ((dif1.X != 0) && (dif2.X != 0))
            {
                m1 = (double)(dif1.Y) / (dif1.X);
                m2 = (double)(dif2.Y) / (dif2.X);
                if (m1 != m2)
                {
                    return false;
                }
                return true;
            }
            return false;
        }

        public bool nearlyParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            return pNearlyParallel(p1, p2, p3, p4);
        }

        bool pNearlyParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
        {
            double m1, m2;
            GeoLibPointF dif1, dif2;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            dif2 = new GeoLibPointF(p4.X - p3.X, p4.Y - p3.Y);
            if ((dif1.Y != 0) && (dif2.Y != 0))
            {
                m1 = dif1.X / dif1.Y;
                m2 = dif2.X / dif2.Y;
                m2 = m2 - m1;
                if ((m2 < 0) && (-0.005 > m2)) { return false; }
                if ((m2 > 0) && (0.005 < m2)) { return false; }
                return true;
            }
            else if ((dif1.X != 0) && (dif2.X != 0))
            {
                m1 = dif1.Y / dif1.X;
                m2 = dif2.Y / dif2.X;
                m2 = m2 - m1;
                if ((m2 < 0) && (-0.005 > m2)) { return false; }
                if ((m2 > 0) && (0.005 < m2)) { return false; }
                return true;
            }
            return false;
        }

        public bool onLine2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            return pOnLine2(p1, p2, p3);
        }

        bool pOnLine2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
        {
            double m1, b1;
            GeoLibPointF dif1;
            int help;
            dif1 = new GeoLibPointF(p2.X - p1.X, p2.Y - p1.Y);
            if ((dif1.X != 0))
            {
                m1 = dif1.Y / dif1.X;
                b1 = -m1 * p1.X + p1.Y;
                if (p3.Y != m1 * p3.X + b1) { return false; }
            }
            else
            {
                if (p3.X != p1.X) { return false; }
            }
            if (p3 == p1) { return false; }
            if (p3 == p2) { return false; }
            if (p1.X > p2.X) { help = p1.X; p1.X = (p2.X); p2.X = (help); }
            if (p1.Y > p2.Y) { help = p1.Y; p1.Y = (p2.Y); p2.Y = (help); }
            if (p3.X < p1.X) { return false; }
            if (p3.Y < p1.Y) { return false; }
            if (p3.X > p2.X) { return false; }
            if (p3.Y > p2.Y) { return false; }

            return true;
        }

        public static int round(double d)
        {
            if (d > 0)
            {
                return (int)Math.Floor(d + 0.5);
            }
            else
            {
                return (int)Math.Floor(d - 0.5);
            }
        }

        public GeoLibPoint round(GeoLibPoint point, int i)
        {
            return pRound(point, i);
        }

        GeoLibPoint pRound(GeoLibPoint point, int i)
        {
            GeoLibPoint p = new GeoLibPoint();
            if (i > 1)
            {
                if (point.X > 0)
                {
                    p.X = ((point.X + (i / 2)) / i) * i;
                }
                else
                {
                    p.X = (((point.X - (i / 2)) / i) * i);
                }
                if (point.Y > 0)
                {
                    p.Y = (((point.Y + (i / 2)) / i) * i);
                }
                else
                {
                    p.Y = (((point.Y - (i / 2)) / i) * i);
                }
            }
            return p;
        }

        public virtual GCCell depend()
        {
            return null;
        }

        public GeoLibPoint[] ellipse(GeoLibPoint center, double radius, double anglestep)
        {
            return pEllipse(center, radius, anglestep);
        }

        GeoLibPoint[] pEllipse(GeoLibPoint center, double radius, double anglestep)
        {
            List<GeoLibPoint> array = new List<GeoLibPoint>();
            double radius1 = radius;
            double radius2 = radius;
            for (double ang = 0; ang < 360; ang += anglestep)
            {
                array.Add(new GeoLibPoint(center.X + (radius1 * Math.Cos(ang / 180 * Math.PI)), center.Y + (radius2 * Math.Sin(ang / 180 * Math.PI))));
            }
            return array.ToArray();
        }

        public virtual void minimum(GeoLibPoint p) { }
        public virtual void maximum(GeoLibPoint p) { }
        public virtual void deleteSelect() { }
        public virtual void moveSelect(GeoLibPoint p) { }
        public virtual void move(GeoLibPoint p) { }
        public virtual void moveToLayer(Int32 layer) { layer_nr = layer; }
        public virtual void moveToDataType(Int32 datatype) { datatype_nr = datatype; }
        public virtual void moveToLayerSelect(Int32 layer) { }
        public virtual void moveToDataTypeSelect(Int32 datatype) { }
        public virtual void mapSelect(GCStrans s) { }
        public virtual void map(GCStrans s) { }
        public virtual void resize(double factor) { }
        public virtual void clean() { }
        public virtual bool correct() { return true; }
        public virtual List<GCElement> cutSelect(GeoLibPoint a, GeoLibPoint b) { return null; }

        // Box
        public virtual bool isBox() { return false; }
        public virtual GCBox getBox() { return null; }
        public virtual GCBox convertToBox() { return null; }

        // Polygon
        public virtual bool isPolygon() { return false; }
        public virtual GCPolygon getPolygon() { return null; }
        public virtual GCPolygon convertToPolygon() { return null; }
        public virtual bool mergeSelect(GCPolygon p) { return false; }
        public virtual void edgeRemoveSelect(int edge) { }
        public virtual bool add(GeoLibPoint[] poly) { return false; }

        // CellRef
        public virtual bool isCellref() { return false; }

        // CellRefArray
        public virtual bool isCellrefArray() { return false; }
        public virtual GCCell getCellref() { return null; }
        public virtual GCCellRefArray getCellrefArray() { return null; }
        public virtual void setMirrorx() { }
        public virtual void clearMirrorx() { }
        public virtual void toggleMirrorx() { }
        public virtual void rotate(double angle) { }
        public virtual void rotate(double angle, GeoLibPoint pos) { }
        public virtual void scale(double factor) { }
        public virtual void scale(GeoLibPoint origin, double factor) { }
        public virtual void setCellRef(GCCell cellRef) { }
        public virtual void setPos(GeoLibPoint p) { }

        public virtual GeoLibPoint getPos() { return null; }
        public virtual List<GCElement> flatSelect() { return null; }
        // and txt;
        public virtual void setName(String s) { }
        public virtual String getName() { return null; }

        // Text
        public virtual bool isText() { return false; }
        public virtual GCTxt getText() { return null; }
        public virtual void setPresentation(int presentation) { }
        public virtual List<GCElement> convertText() { return null; }

        // Path
        public virtual bool isPath() { return false; }
        public virtual GCPath getPath() { return null; }
        public virtual void expandCaps(double a, double b) { }

        // Path and text
        public virtual void setWidth(int width) { }
        public virtual int getWidth() { return 0; }
        public virtual GCPolygon closeToPolygon() { return null; }
        public virtual bool mergeSelect(GCPath aPath) { return false; }
        public virtual void setCap(int v) { }
        public virtual void nextPoint() { }

        // IO
        public virtual void saveGDS(gdsWriter gw) { }
        public virtual void saveOASIS(oasWriter ow) { }
    }
}
