using gds;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;
using geoWrangler;

namespace geoCoreLib;

public class GCCell
{
    public string cellName { get; set; }

    public bool saved { get; set; }

    public short modyear { get; set; }
    public short modmonth { get; set; }
    public short modday { get; set; }
    public short modhour { get; set; }
    public short modmin { get; set; }
    public short modsec { get; set; }

    public short accyear { get; set; }
    public short accmonth { get; set; }
    public short accday { get; set; }
    public short acchour { get; set; }
    public short accmin { get; set; }
    public short accsec { get; set; }

    public List<GCElement> elementList { get; set; }

    public GCCell()
    {
        pGCCell();
    }

    private void pGCCell()
    {
        DateTime now = DateTime.Now;
        modyear = (short)now.Year;
        modmonth = (short)now.Month;
        modday = (short)now.Day;
        modhour = (short)now.Hour;
        modmin = (short)now.Minute;
        modsec = (short)now.Second;
        accyear = (short)now.Year;
        accmonth = (short)now.Month;
        accday = (short)now.Day;
        acchour = (short)now.Hour;
        accmin = (short)now.Minute;
        accsec = (short)now.Second;
        cellName = "noname";
        elementList = new List<GCElement>();
    }
    
    public void setName(string name)
    {
        cellName = name;
    }
    
    public string getName()
    {
        return cellName;
    }

    public void addBox(int x, int y, int b, int h, int layer, int datatype)
    {
        pAddBox(x, y, b, h, layer, datatype);
    }

    private void pAddBox(int x, int y, int b, int h, int layer, int datatype)
    {
        addElement(new GCBox(x, y, b, h, layer, datatype));
    }

    public void addBox(Path64 points, int layer, int datatype)
    {
        pAddBox(points, layer, datatype);
    }

    private void pAddBox(Path64 points, int layer, int datatype)
    {
        GCElement e;

        if (points.Count != 5)
        {
            e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
            return;
        }
        if ((points[0].X != points[4].X) || (points[0].Y != points[4].Y))
        {
            e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
            return;
        }

        int x2, y2;
        Point64 p = new(points[0]);
        int x1 = (int)p.X;
        int y1 = (int)p.Y;
        p = new(points[1]);
        if (p.X < x1)
        {
            x2 = x1;
            x1 = (int)p.X;
        }
        else
        {
            x2 = (int)p.X;
        }
        if (p.Y < y1)
        {
            y2 = y1;
            y1 = (int)p.Y;
        }
        else
        {
            y2 = (int)p.Y;
        }
        for (int i = 2; i < 4; i++)
        {
            p = new(points[i]);
            if (p.X < x1)
            {
                x1 = (int)p.X;
            }
            if (p.X > x2)
            {
                x2 = (int)p.X;
            }
            if (p.Y < y1)
            {
                y1 = (int)p.Y;
            }
            if (p.Y > y2)
            {
                y2 = (int)p.Y;
            }
        }
        bool b = true;
        for (int i = 0; i < 4; i++)
        {
            p = new(points[i]);
            if (p.X != x1 && p.X != x2)
            {
                b = false;
            }

            if (p.Y != y1 && p.Y != y2)
            {
                b = false;
            }
        }

        switch (b)
        {
            case true:
                int box_width = x2 - x1;
                int box_height = y2 - y1;
                int box_x = x1;
                int box_y = y1;
                e = new GCBox(box_x, box_y, box_width, box_height, layer, datatype);
                pAddElement(e);
                break;
            default:
                e = new GCPolygon(points, layer, datatype);
                pAddElement(e);
                break;
        }
    }

    public void addCircle(int layer, int datatype, Point64 center, double radius)
    {
        pAddCircle(layer, datatype, center, radius);
    }

    private void pAddCircle(int layer, int datatype, Point64 center, double radius)
    {
        GCElement e = new();
        Path64 points = e.ellipse(center, radius, GCSetup.circularDefault);
        e = new GCPolygon(points, layer, datatype);
        pAddElement(e);
    }

    public void addElement(GCElement element)
    {
        pAddElement(element);
    }

    private void pAddElement(GCElement element)
    {
        elementList.Add(element);
    }

    public void addPath(Path64 points, int layer, int datatype)
    {
        pAddPath(points, layer, datatype);
    }

    private void pAddPath(Path64 points, int layer, int datatype)
    {
        GCElement e = new GCPath(points, layer, datatype);
        pAddElement(e);
    }

    public void addPolygon(Path64 points, int layer, int datatype)
    {
        pAddPolygon(points, layer, datatype);
    }

    private void pAddPolygon(Path64 points, int layer, int datatype)
    {
        // Might this be a box?
        bool box = false;
        // Box is defined by 5 points
        if (points.Count == 5)
        {
            // Check that the distances are equi-distant on each side, using the diagonal opposite
            if (Math.Abs(points[0].X - points[2].X) == Math.Abs(points[0].Y - points[2].Y))
            {
                double[] angles = GeoWrangler.angles(points, true).Distinct().ToArray();
                if (angles.Length == 1)
                {
                    box = true;
                }
            }
        }

        GCElement e;
        if (box)
        {
            pAddBox(points, layer, datatype);
        }
        else
        {
            e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
        }
    }

    public void addCellref(GCCell c, Point64 pos)
    {
        pAddCellref(c, pos, 0, 1, false);
    }

    public void addCellref(GCCell c, Point64 pos, double angle, double mag, bool mirror_x)
    {
        pAddCellref(c, pos, angle, mag, mirror_x);
    }

    private void pAddCellref(GCCell c, Point64 pos, double angle, double mag, bool mirror_x)
    {
        GCElement e = new GCCellref(c, pos, angle, mag, mirror_x);
        pAddElement(e);
    }
    
    public void addCellrefArray(GCCell c, Point64 pos, Repetition r)
    {
        pAddCellrefArray(c, pos, 0, 1, false, r);
    }

    public void addCellrefArray(GCCell c, Point64 pos, double angle, double mag, bool mirror_x, Repetition r)
    {
        pAddCellrefArray(c, pos, angle, mag, mirror_x, r);
    }

    private void pAddCellrefArray(GCCell c, Point64 pos, double angle, double mag, bool mirror_x, Repetition r)
    {
        GCElement e = new GCCellRefArray(c, angle, mag, mirror_x, r);
        e.move(pos);
        pAddElement(e);
    }

    public void addCellrefArray(GCCell c, Path64 array, int anzx, int anzy)
    {
        pAddCellrefArray(c, array, anzx, anzy);
    }

    private void pAddCellrefArray(GCCell c, Path64 array, int anzx, int anzy)
    {
        GCElement e = new GCCellRefArray(c, array, anzx, anzy);
        pAddElement(e);
    }

    public void addCellrefArray(GCCell c, Point64 pos1, Point64 pos2, int anzx, int anzy)
    {
        pAddCellrefArray(c, pos1, pos2, anzx, anzy);
    }

    private void pAddCellrefArray(GCCell c, Point64 pos1, Point64 pos2, int anzx, int anzy)
    {
        GCElement e = new GCCellRefArray(c, pos1, pos2, anzx, anzy);
        pAddElement(e);
    }

    public void addCellrefs(int count)
    {
        pAddCellrefs(count);
    }

    private void pAddCellrefs(int count)
    {
        for (int i = 0; i < count; i++)
        {
            elementList.Add(new GCCellref());
        }
    }

    public void addCellref()
    {
        pAddCellref();
    }

    private void pAddCellref()
    {
        GCElement e = new GCCellref();
        pAddElement(e);
    }

    public void addText(int layer, int datatype, Point64 pos, string text)
    {
        pAddText(layer, datatype, pos, text);
    }

    private void pAddText(int layer, int datatype, Point64 pos, string text)
    {
        GCElement e = new GCTxt(layer, datatype, pos, text);
        pAddElement(e);
    }

    public bool depend(GCCell cell)
    {
        return pDepend(cell);
    }

    private bool pDepend(GCCell cell)
    {
        bool b = false;
        foreach (GCCell cellhelp in from t in elementList where t != null select t.depend())
        {
            if (cellhelp == cell)
            {
                b = true;
            }
            else
            {
                if (cellhelp != null && !b)
                {
                    b = cellhelp.depend(cell);
                }
            }
        }
        return b;
    }

    public void move(Point64 p)
    {
        pMove(p);
    }

    private void pMove(Point64 p)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].move(p);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void moveSelect(Point64 p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(Point64 p)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveSelect(p);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void moveToDataType(int datatype)
    {
        pMoveToDataType(datatype);
    }

    private void pMoveToDataType(int datatype)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToDataType(datatype);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void moveToLayer(int layer)
    {
        pMoveToLayer(layer);
    }

    private void pMoveToLayer(int layer)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToLayer(layer);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void moveToDataTypeSelect(int datatype)
    {
        pMoveToDataTypeSelect(datatype);
    }

    private void pMoveToDataTypeSelect(int datatype)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToDataTypeSelect(datatype);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void moveToLayerSelect(int layer)
    {
        pMoveToLayerSelect(layer);
    }

    private void pMoveToLayerSelect(int layer)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToLayerSelect(layer);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void resize(double size)
    {
        pResize(size);
    }

    private void pResize(double size)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    elementList[f].resize(size);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void rotateSelect(double angle, Point64 pos)
    {
        pRotateSelect(angle, pos);
    }

    private void pRotateSelect(double angle, Point64 pos)
    {
        GCStrans m = new();
        m.translate(pos.X, pos.Y);
        m.rotate(angle);
        m.translate(-pos.X, -pos.Y);
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] == null)
                {
                    return;
                }

                if (elementList[f].isBox() && elementList[f].select)
                {
                    GCPolygon p = elementList[f].convertToPolygons(1.0)[0];
                    p.map(m);
                    GCBox b = p.convertToBox();
                    if (b != null)
                    {
                        elementList[f] = b;
                    }
                    else
                    {
                        elementList[f] = p;
                    }
                    elementList[f].select = true;
                }
                else
                {
                    elementList[f].mapSelect(m);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void scaleSelect(Point64 pos, Point64 p2, Point64 p3)
    {
        pScaleSelect(pos, p2, p3);
    }

    private void pScaleSelect(Point64 pos, Point64 p2, Point64 p3)
    {
        GCStrans m = new();
        m.translate(pos.X, pos.Y);
        double x = (double)(p3.X - pos.X) / (p2.X - pos.X);
        double y = (double)(p3.Y - pos.Y) / (p2.Y - pos.Y);
        x = x switch
        {
            0 => 1,
            _ => x
        };
        y = y switch
        {
            0 => 1,
            _ => y
        };
        x = (p2.X - pos.X) switch
        {
            0 => 1,
            _ => x
        };
        y = (p2.Y - pos.Y) switch
        {
            0 => 1,
            _ => y
        };
        m.scale(x, y);
        m.translate(-pos.X, -pos.Y);
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] == null)
                {
                    return;
                }

                if (elementList[f].isBox() && elementList[f].select)
                {
                    GCPolygon p = elementList[f].convertToPolygons(1.0)[0];
                    p.map(m);
                    GCBox b = p.convertToBox();
                    if (b != null)
                    {
                        elementList[f] = b;
                    }
                    else
                    {
                        elementList[f] = p;
                    }
                    elementList[f].select = true;
                }
                else
                {
                    elementList[f].mapSelect(m);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void mirrorSelect(Point64 p1, Point64 p2)
    {
        pMirrorSelect(p1, p2);
    }

    private void pMirrorSelect(Point64 p1, Point64 p2)
    {
        GCStrans m = new();
        if (p1 == p2)
        {
            m.translate(p1.X, p1.Y);
            m.scale(-1);
            m.translate(-p1.X, -p1.Y);
        }
        else
        {
            double angle;
            if ((p1.X - p2.X) == 0)
            {
                angle = 90;
            }
            else
            {
                angle = Math.Atan((p1.Y - p2.Y) / (double)(p1.X - p2.X)) / 2 / Math.PI * 360;
            }

            m.translate(p1.X, p1.Y);
            m.rotate(angle);
            m.toggleMirror_x();
            m.rotate(-angle);
            m.translate(-p1.X, -p1.Y);
        }
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] == null)
                {
                    return;
                }

                if (elementList[f].isBox() && elementList[f].select)
                {
                    GCPolygon p = elementList[f].convertToPolygons(1.0)[0];
                    p.map(m);
                    GCBox b = p.convertToBox();
                    if (b != null)
                    {
                        elementList[f] = b;
                    }
                    else
                    {
                        elementList[f] = p;
                    }
                    elementList[f].select = true;
                }
                else
                {
                    elementList[f].mapSelect(m);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }
        
    public void minimum(ref Point64 pos)
    {
        pMinimum(ref pos);
    }

    private void pMinimum(ref Point64 pos)
    {
        foreach (GCElement el in elementList)
        {
            if (el != null)
            {
                el.minimum(ref pos);
            }
        }
    }

    public void maximum(ref Point64 pos)
    {
        pMaximum(ref pos);
    }

    private void pMaximum(ref Point64 pos)
    {
        foreach (GCElement el in elementList)
        {
            el?.maximum(ref pos);
        }
    }

    public void clean()
    {
        pClean();
    }

    private void pClean()
    {
        for (int f = 0; f < elementList.Count; f++)
        {
            if (elementList[f] == null)
            {
                continue;
            }

            if (!elementList[f].correct())
            {
                elementList[f] = null;
            }
        }
        elementList.RemoveAll(null!);
    }

    public void updateCellref(GCCell oldCell, GCCell newCell)
    {
        pUpdateCellref(oldCell, newCell);
    }

    private void pUpdateCellref(GCCell oldCell, GCCell newCell)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, f =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] == null)
                {
                    return;
                }

                if (elementList[f].depend() == oldCell)
                {
                    elementList[f].setCellRef(newCell);
                }
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public bool dependNotSaved()
    {
        return pDependNotSaved();
    }

    private bool pDependNotSaved()
    {
        switch (saved)
        {
            case true:
                return false;
        }

        return elementList.Select(element => element?.depend()).Any(cellHelp => cellHelp is {saved: false});
    }

    public void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        //bgnstr
        gw.bw.Write((ushort)28);
        gw.bw.Write(gdsValues.sBGNSTR);
        // Get date and time.

        // Modification
        gw.bw.Write(modyear);
        gw.bw.Write(modmonth);
        gw.bw.Write(modday);
        gw.bw.Write(modhour);
        gw.bw.Write(modmin);
        gw.bw.Write(modsec);

        // Access
        gw.bw.Write(accyear);
        gw.bw.Write(accmonth);
        gw.bw.Write(accday);
        gw.bw.Write(acchour);
        gw.bw.Write(accmin);
        gw.bw.Write(accsec);
        gw.writeString(cellName, 6);
        foreach (GCElement el in elementList)
        {
            if (el != null)
            {
                el.saveGDS(gw);
            }
        }
        //endstr
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sENDSTR);
        saved = true;
    }

    public void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    private void pSaveOASIS(oasWriter ow)
    {
        foreach (GCElement el in elementList)
        {
            if (el != null)
            {
                el.saveOASIS(ow);
            }
        }

        saved = true;
    }

    public List<GCPolygon> convertToPolygons(double scaleFactor, int layer = -1, int datatype = -1)
    {
        return pConvertToPolygons(scaleFactor, layer, datatype);
    }

    private List<GCPolygon> pConvertToPolygons(double scaleFactor, int layer = -1, int datatype = -1)
    {
        List<GCPolygon> ret = new();
        foreach (GCElement t in elementList.Where(t => t != null))
        {
            if (layer == -1 || datatype == -1)
            {
                ret.AddRange(t.convertToPolygons(scaleFactor));
            }
            else
            {
                ret.AddRange(t.convertToPolygons(scaleFactor).Where(p => p.layer_nr == layer && p.datatype_nr == datatype));
            }
        }

        return ret;
    }
}