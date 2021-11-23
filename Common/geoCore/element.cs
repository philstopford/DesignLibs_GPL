using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;

namespace geoCoreLib;

public class GCElement
{
    public int layer_nr { get; set; }
    public int datatype_nr { get; set; }
    public bool select { get; set; }

    public double angle(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        return pAngle(p1, p2, p3);
    }

    private double pAngle(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        double a1, a2;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        GeoLibPointF dif2 = new(p3.X - p2.X, p3.Y - p2.Y);
        if (dif1.X != 0)
        {
            a1 = Math.Atan(dif1.Y / dif1.X) / 2 / Math.PI * 360;
            switch (dif1.X)
            {
                case < 0:
                    a1 -= 180;
                    break;
            }

            switch (a1)
            {
                case < -180:
                    a1 += 360;
                    break;
            }
        }
        else
        {
            a1 = dif1.Y switch
            {
                > 0 => 90,
                _ => -90
            };
        }
        if (dif2.X != 0)
        {
            a2 = Math.Atan(dif2.Y / dif2.X) / 2 / Math.PI * 360;
            switch (dif2.X)
            {
                case < 0:
                    a2 -= 180;
                    break;
            }

            switch (a2)
            {
                case < -180:
                    a2 += 360;
                    break;
            }
        }
        else
        {
            a2 = dif2.Y switch
            {
                > 0 => 90,
                _ => -90
            };
        }
        a1 = a2 - a1;
        switch (a1)
        {
            case <= -180:
                a1 += 360;
                break;
        }

        switch (a1)
        {
            case >= 180:
                a1 -= 360;
                break;
        }
        return a1;
    }

    public double angle(GeoLibPoint p1, GeoLibPoint p2)
    {
        return pAngle(p1, p2);
    }

    private double pAngle(GeoLibPoint p1, GeoLibPoint p2)
    {
        double a1;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        if (dif1.X != 0)
        {
            a1 = Math.Atan(dif1.Y / dif1.X) / 2 / Math.PI * 360;
            switch (dif1.X)
            {
                case < 0:
                    a1 -= 180;
                    break;
            }

            switch (a1)
            {
                case < -180:
                    a1 += 360;
                    break;
            }
        }
        else
        {
            a1 = dif1.Y switch
            {
                > 0 => 90,
                _ => -90
            };
        }

        switch (a1)
        {
            case <= -180:
                a1 += 360;
                break;
        }

        switch (a1)
        {
            case > 180:
                a1 -= 360;
                break;
        }
        return a1;
    }

    public bool cutPoint2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4, GeoLibPoint pc)
    {
        return pCutPoint2(p1, p2, p3, p4, pc);
    }

    private bool pCutPoint2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4, GeoLibPoint pc)
    {
        //between  p1-p2+p3-p4
        double m1, m2, b1, b2;
        int help;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        GeoLibPointF dif2 = new(p4.X - p3.X, p4.Y - p3.Y);
        pc.X = 1 << 30;
        pc.Y = 1 << 30;
        if (dif1.X != 0 && dif2.X != 0)
        {
            m1 = dif1.Y / dif1.X;
            m2 = dif2.Y / dif2.X;
            b1 = -m1 * p1.X + p1.Y;
            b2 = -m2 * p3.X + p3.Y;
            switch (Math.Abs(m1 - m2))
            {
                case <= double.Epsilon:
                    return false;
            }
            pc.X = round((b2 - b1) / (m1 - m2));
            switch (dif1.Y)
            {
                case 0:
                    pc.Y = p1.Y;
                    break;
                default:
                {
                    pc.Y = dif2.Y switch
                    {
                        0 => p3.Y,
                        _ => round(m1 * pc.X + b1)
                    };

                    break;
                }
            }
        }
        switch (dif1.X)
        {
            case 0 when dif2.X != 0:
            {
                m2 = dif2.Y / dif2.X;
                b2 = -m2 * p3.X + p3.Y;
                pc.X = p1.X;
                pc.Y = dif2.Y != 0 ? round(m2 * pc.X + b2) : p3.Y;

                break;
            }
        }
        if (dif1.X != 0 && dif2.X == 0)
        {
            m1 = dif1.Y / dif1.X;
            b1 = -m1 * p1.X + p1.Y;
            pc.X = p3.X;
            pc.Y = dif2.Y != 0 ? round(m1 * pc.X + b1) : p1.Y;
        }

        switch (dif1.X)
        {
            case 0 when dif2.X == 0:
                return false;
        }

        if (pc == p1 || pc == p2 || pc == p3 || pc == p4)
        {
            return false;
        }

        if (p1.X > p2.X)
        {
            help = p1.X;
            p1.X = p2.X;
            p2.X = help;
        }

        if (p3.X > p4.X)
        {
            help = p3.X;
            p3.X = p4.X;
            p4.X = help;
        }

        if (pc.X < p1.X || pc.X > p2.X || pc.X < p3.X || pc.X > p4.X)
        {
            return false;
        }

        if (p1.Y > p2.Y)
        {
            help = p1.Y;
            p1.Y = p2.Y;
            p2.Y = help;
        }

        if (p3.Y > p4.Y)
        {
            help = p3.Y;
            p3.Y = p4.Y;
            p4.Y = help;
        }

        return pc.Y >= p1.Y && pc.Y <= p2.Y && pc.Y >= p3.Y && pc.Y <= p4.Y;
    }

    public double distance(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        return pDistance(p1, p2, p3);
    }

    private double pDistance(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        // >0 if left off line
        GeoLibPoint dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        if (dif1.X == 0)
        {
            return dif1.Y switch
            {
                > 0 => p1.X - p3.X,
                _ => p3.X - p1.X
            };
        }

        double m1 = (double)dif1.Y / dif1.X;
        double b1 = -m1 * p1.X + p1.Y;
        double wert = m1 * p3.X + b1 - p3.Y;
        return dif1.Y switch
        {
            0 when dif1.X > 0 => -wert,
            0 => wert,
            _ => dif1.X switch
            {
                > 0 => -wert,
                < 0 => +wert,
                _ => wert
            }
        };

    }

    public double distance(GeoLibPoint p1, GeoLibPoint p2)
    {
        return pDistance(p1, p2);
    }

    private double pDistance(GeoLibPoint p1, GeoLibPoint p2)
    {
        int dx = p1.X - p2.X;
        int dy = p1.Y - p2.Y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    public bool identical(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        return pIdentical(p1, p2, p3, p4);
    }

    private bool pIdentical(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        double m1, m2, b1, b2;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        GeoLibPointF dif2 = new(p4.X - p3.X, p4.Y - p3.Y);
        if (dif1.Y != 0 && dif2.Y != 0)
        {
            m1 = dif1.X / dif1.Y;
            m2 = dif2.X / dif2.Y;

            switch (Math.Abs(m1 - m2))
            {
                case > double.Epsilon:
                    return false;
            }
            b1 = -m1 * p1.Y + p1.X;
            b2 = -m2 * p3.Y + p3.X;

            return Math.Abs(b1 - b2) switch
            {
                > double.Epsilon => false,
                _ => true
            };
        }

        if (dif1.X == 0 || dif2.X == 0)
        {
            return false;
        }

        m1 = dif1.Y / dif1.X;
        m2 = dif2.Y / dif2.X;

        switch (Math.Abs(m1 - m2))
        {
            case > double.Epsilon:
                return false;
        }
        b1 = -m1 * p1.X + p1.Y;
        b2 = -m2 * p3.X + p3.Y;

        return Math.Abs(b1 - b2) switch
        {
            > double.Epsilon => false,
            _ => true
        };
    }

    public double length(GeoLibPoint p)
    {
        return pLength(p);
    }

    private double pLength(GeoLibPoint p)
    {
        return Math.Sqrt(p.X * p.X + p.Y * p.Y);
    }

    public bool parallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        return pParallel(p1, p2, p3, p4);
    }

    private bool pParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        double m1, m2;
        GeoLibPoint dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        GeoLibPoint dif2 = new(p4.X - p3.X, p4.Y - p3.Y);
        if (dif1.Y != 0 && dif2.Y != 0)
        {
            m1 = (double)dif1.X / dif1.Y;
            m2 = (double)dif2.X / dif2.Y;
            return Math.Abs(m1 - m2) switch
            {
                > double.Epsilon => false,
                _ => true
            };
        }

        if (dif1.X == 0 || dif2.X == 0)
        {
            return false;
        }

        m1 = (double)dif1.Y / dif1.X;
        m2 = (double)dif2.Y / dif2.X;
        return Math.Abs(m1 - m2) switch
        {
            > double.Epsilon => false,
            _ => true
        };
    }

    public bool nearlyParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        return pNearlyParallel(p1, p2, p3, p4);
    }

    private bool pNearlyParallel(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3, GeoLibPoint p4)
    {
        double m1, m2;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        GeoLibPointF dif2 = new(p4.X - p3.X, p4.Y - p3.Y);
        if (dif1.Y != 0 && dif2.Y != 0)
        {
            m1 = dif1.X / dif1.Y;
            m2 = dif2.X / dif2.Y;
            m2 -= m1;
            switch (m2)
            {
                case < 0 and < -0.005:
                case > 0 and > 0.005:
                    return false;
                default:
                    return true;
            }
        }

        if (dif1.X == 0 || dif2.X == 0)
        {
            return false;
        }

        m1 = dif1.Y / dif1.X;
        m2 = dif2.Y / dif2.X;
        m2 -= m1;
        switch (m2)
        {
            case < 0 and < -0.005:
            case > 0 and > 0.005:
                return false;
            default:
                return true;
        }
    }

    public bool onLine2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        return pOnLine2(p1, p2, p3);
    }

    private bool pOnLine2(GeoLibPoint p1, GeoLibPoint p2, GeoLibPoint p3)
    {
        int help;
        GeoLibPointF dif1 = new(p2.X - p1.X, p2.Y - p1.Y);
        if (dif1.X != 0)
        {
            double m1 = dif1.Y / dif1.X;
            double b1 = -m1 * p1.X + p1.Y;
            switch (Math.Abs(p3.Y - (m1 * p3.X + b1)))
            {
                case > double.Epsilon:
                    return false;
            }
        }
        else
        {
            if (p3.X != p1.X) { return false; }
        }

        if (p3 == p1 || p3 == p2)
        {
            return false;
        }

        if (p1.X > p2.X)
        {
            help = p1.X;
            p1.X = p2.X;
            p2.X = help;
        }

        if (p1.Y > p2.Y)
        {
            help = p1.Y;
            p1.Y = p2.Y;
            p2.Y = help;
        }

        return p3.X >= p1.X && p3.Y >= p1.Y && p3.X <= p2.X && p3.Y <= p2.Y;
    }

    public static int round(double d)
    {
        return d switch
        {
            > 0 => (int) Math.Floor(d + 0.5),
            _ => (int) Math.Floor(d - 0.5)
        };
    }

    public GeoLibPoint round(GeoLibPoint point, int i)
    {
        return pRound(point, i);
    }

    private GeoLibPoint pRound(GeoLibPoint point, int i)
    {
        GeoLibPoint p = new();
        switch (i)
        {
            case > 1:
            {
                p.X = point.X switch
                {
                    > 0 => (point.X + i / 2) / i * i,
                    _ => (point.X - i / 2) / i * i
                };
                p.Y = point.Y switch
                {
                    > 0 => (point.Y + i / 2) / i * i,
                    _ => (point.Y - i / 2) / i * i
                };

                break;
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

    private GeoLibPoint[] pEllipse(GeoLibPoint center, double radius, double anglestep)
    {
        List<GeoLibPoint> array = new();
        double radius1 = radius;
        double radius2 = radius;
        for (double ang = 0; ang < 360; ang += anglestep)
        {
            array.Add(new GeoLibPoint(center.X + radius1 * Math.Cos(ang / 180 * Math.PI), center.Y + radius2 * Math.Sin(ang / 180 * Math.PI)));
        }
        return array.ToArray();
    }

    public virtual void minimum(GeoLibPoint p) { }
    public virtual void maximum(GeoLibPoint p) { }
    public virtual void deleteSelect() { }
    public virtual void moveSelect(GeoLibPoint p) { }
    public virtual void move(GeoLibPoint p) { }
    public virtual void moveToLayer(int layer) { layer_nr = layer; }
    public virtual void moveToDataType(int datatype) { datatype_nr = datatype; }
    public virtual void moveToLayerSelect(int layer) { }
    public virtual void moveToDataTypeSelect(int datatype) { }
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
    public virtual List<GCPolygon> convertToPolygons() { return null; }
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
    public virtual void setName(string s) { }
    public virtual string getName() { return null; }

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