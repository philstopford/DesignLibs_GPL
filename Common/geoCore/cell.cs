using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

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
        modyear = 0;
        modmonth = 0;
        modday = 0;
        modhour = 0;
        modmin = 0;
        modsec = 0;
        accyear = 0;
        accmonth = 0;
        accday = 0;
        acchour = 0;
        accmin = 0;
        accsec = 0;
        cellName = "noname";
        elementList = new List<GCElement>();
    }

    public void addBox(int x, int y, int b, int h, int layer, int datatype)
    {
        pAddBox(x, y, b, h, layer, datatype);
    }

    private void pAddBox(int x, int y, int b, int h, int layer, int datatype)
    {
        addElement(new GCBox(x, y, b, h, layer, datatype));
    }

    public void addBox(GeoLibPoint[] points, int layer, int datatype)
    {
        pAddBox(points, layer, datatype);
    }

    private void pAddBox(GeoLibPoint[] points, int layer, int datatype)
    {
        GCElement e;

        if (points.Length != 5)
        {
            e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
            return;
        }
        if (points[0] != points[4])
        {
            e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
            return;
        }
        GeoLibPoint p;
        int x2, y2;
        p = points[0];
        int x1 = p.X;
        int y1 = p.Y;
        p = points[1];
        if (p.X < x1)
        {
            x2 = x1;
            x1 = p.X;
        }
        else
        {
            x2 = p.X;
        }
        if (p.Y < y1)
        {
            y2 = y1;
            y1 = p.Y;
        }
        else
        {
            y2 = p.Y;
        }
        for (int i = 2; i < 4; i++)
        {
            p = points[i];
            if (p.X < x1)
            {
                x1 = p.X;
            }
            if (p.X > x2)
            {
                x2 = p.X;
            }
            if (p.Y < y1)
            {
                y1 = p.Y;
            }
            if (p.Y > y2)
            {
                y2 = p.Y;
            }
        }
        bool b = true;
        for (int i = 0; i < 4; i++)
        {
            p = points[i];
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
                e = new GCBox(x1, y1, x2 - x1 + 1, y2 - y1 + 1, layer, datatype);
                pAddElement(e);
                break;
            default:
                e = new GCPolygon(points, layer, datatype);
                pAddElement(e);
                break;
        }
    }

    public GCElement addCircle(int layer, int datatype, GeoLibPoint center, double radius)
    {
        return pAddCircle(layer, datatype, center, radius);
    }

    private GCElement pAddCircle(int layer, int datatype, GeoLibPoint center, double radius)
    {
        GCElement e = new();
        GeoLibPoint[] points = e.ellipse(center, radius, GCSetup.circularDefault);
        e = new GCPolygon(points, layer, datatype);
        return e;
    }

    public void addElement(GCElement element)
    {
        pAddElement(element);
    }

    private void pAddElement(GCElement element)
    {
        elementList.Add(element);
    }

    public void addPath(GeoLibPoint[] points, int layer, int datatype)
    {
        pAddPath(points, layer, datatype);
    }

    private void pAddPath(GeoLibPoint[] points, int layer, int datatype)
    {
        GCElement e = new GCPath(points, layer, datatype);
        pAddElement(e);
    }

    public void addPolygon(GeoLibPoint[] points, int layer, int datatype)
    {
        pAddPolygon(points, layer, datatype);
    }

    private void pAddPolygon(GeoLibPoint[] points, int layer, int datatype)
    {
        GCElement e = new GCPolygon(points, layer, datatype);
        pAddElement(e);
    }

    public void addCellref(GCCell c, GeoLibPoint pos)
    {
        pAddCellref(c, pos);
    }

    private void pAddCellref(GCCell c, GeoLibPoint pos)
    {
        GCElement e = new GCCellref(c, pos);
        pAddElement(e);
    }

    public void addCellrefArray(GCCell c, GeoLibPoint[] array, int anzx, int anzy)
    {
        pAddCellrefArray(c, array, anzx, anzy);
    }

    private void pAddCellrefArray(GCCell c, GeoLibPoint[] array, int anzx, int anzy)
    {
        GCElement e = new GCCellRefArray(c, array, anzx, anzy);
        pAddElement(e);
    }

    public void addCellrefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, int anzx, int anzy)
    {
        pAddCellrefArray(c, pos1, pos2, anzx, anzy);
    }

    private void pAddCellrefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, int anzx, int anzy)
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
        GCElement[] e = new GCCellref[count];
        elementList.AddRange(e.ToList());
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

    public void addText(int layer, int datatype, GeoLibPoint pos, string text)
    {
        pAddText(layer, datatype, pos, text);
    }

    private void pAddText(int layer, int datatype, GeoLibPoint pos, string text)
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
        foreach (GCElement t in elementList)
        {
            if (t != null)
            {
                GCCell cellhelp = t.depend();
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
        }
        return b;
    }

    public void move(GeoLibPoint p)
    {
        pMove(p);
    }

    private void pMove(GeoLibPoint p)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, (f) =>
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

    public void moveSelect(GeoLibPoint p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(GeoLibPoint p)
    {
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, (f) =>
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
        Parallel.For(0, elementListCount, (f) =>
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
        Parallel.For(0, elementListCount, (f) =>
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
        Parallel.For(0, elementListCount, (f) =>
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
        Parallel.For(0, elementListCount, (f) =>
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
        Parallel.For(0, elementListCount, (f) =>
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

    public void rotateSelect(double angle, GeoLibPoint pos)
    {
        pRotateSelect(angle, pos);
    }

    private void pRotateSelect(double angle, GeoLibPoint pos)
    {
        GCStrans m = new();
        m.translate(pos.X, pos.Y);
        m.rotate(angle);
        m.translate(-pos.X, -pos.Y);
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, (f) =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].isBox() && elementList[f].@select)
                    {
                        GCPolygon p = elementList[f].convertToPolygons()[0];
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
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void scaleSelect(GeoLibPoint pos, GeoLibPoint p2, GeoLibPoint p3)
    {
        pScaleSelect(pos, p2, p3);
    }

    private void pScaleSelect(GeoLibPoint pos, GeoLibPoint p2, GeoLibPoint p3)
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
        Parallel.For(0, elementListCount, (f) =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].isBox() && elementList[f].@select)
                    {
                        GCPolygon p = elementList[f].convertToPolygons()[0];
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
            }
#if !GCSINGLETHREADED
        );
#endif
    }

    public void mirrorSelect(GeoLibPoint p1, GeoLibPoint p2)
    {
        pMirrorSelect(p1, p2);
    }

    private void pMirrorSelect(GeoLibPoint p1, GeoLibPoint p2)
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
            double angle = (p1.X - p2.X) switch
            {
                0 => 90,
                _ => Math.Atan((p1.Y - p2.Y) / (double) (p1.X - p2.X)) / 2 / Math.PI * 360
            };
            m.translate(p1.X, p1.Y);
            m.rotate(angle);
            m.toggleMirror_x();
            m.rotate(-angle);
            m.translate(-p1.X, -p1.Y);
        }
        int elementListCount = elementList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, elementListCount, (f) =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].isBox() && elementList[f].@select)
                    {
                        GCPolygon p = elementList[f].convertToPolygons()[0];
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
            }
#if !GCSINGLETHREADED
        );
#endif
    }
        
    public void minimum(GeoLibPoint pos)
    {
        pMinimum(pos);
    }

    private void pMinimum(GeoLibPoint pos)
    {
        foreach (GCElement el in elementList)
        {
            if (el != null)
            {
                el.minimum(pos);
            }
        }
    }

    public void maximum(GeoLibPoint pos)
    {
        pMaximum(pos);
    }

    private void pMaximum(GeoLibPoint pos)
    {
        foreach (GCElement el in elementList)
        {
            if (el != null)
            {
                el.maximum(pos);
            }
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
            if (elementList[f] != null)
            {
                if (!elementList[f].correct())
                {
                    elementList[f] = null;
                }
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
        Parallel.For(0, elementListCount, (f) =>
#else
            for (int f = 0; f < elementListCount; f++)
#endif
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].depend() == oldCell)
                    {
                        elementList[f].setCellRef(newCell);
                    }
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
        bool b = false;

        switch (saved)
        {
            case true:
                return false;
        }

        foreach (GCElement element in elementList)
        {
            if (element != null)
            {
                GCCell cellHelp = element.depend();
                if (cellHelp != null && !b)
                {
                    if (!cellHelp.saved)
                    {
                        b = true;
                        break;
                    }
                }
            }
            if (b)
            {
                break;
            }
        }

        return b;
    }

    public void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        //bgnstr
        gw.bw.Write((ushort)28);
        gw.bw.Write((byte)5);
        gw.bw.Write((byte)2);
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
        gw.bw.Write((byte)7);
        gw.bw.Write((byte)0);
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

    public List<GCPolygon> convertToPolygons(int layer = -1, int datatype = -1)
    {
        return pConvertToPolygons(layer, datatype);
    }

    private List<GCPolygon> pConvertToPolygons(int layer = -1, int datatype = -1)
    {
        List<GCPolygon> ret = new();
        for (int f = 0; f < elementList.Count; f++)
        {
            if (elementList[f] != null)
            {
                if (layer == -1 || datatype == -1)
                {
                    ret.AddRange(elementList[f].convertToPolygons());
                }
                else
                {
                    ret.AddRange(elementList[f].convertToPolygons().Where(p => p.layer_nr == layer && p.datatype_nr == datatype));
                }
            }
        }

        return ret;
    }
}