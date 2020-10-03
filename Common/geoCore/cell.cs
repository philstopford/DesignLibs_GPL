using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoCoreLib
{
    public class GCCell
    {
        public string cellName { get; set; }

        public bool saved { get; set; }

        public Int16 modyear { get; set; }
        public Int16 modmonth { get; set; }
        public Int16 modday { get; set; }
        public Int16 modhour { get; set; }
        public Int16 modmin { get; set; }
        public Int16 modsec { get; set; }

        public Int16 accyear { get; set; }
        public Int16 accmonth { get; set; }
        public Int16 accday { get; set; }
        public Int16 acchour { get; set; }
        public Int16 accmin { get; set; }
        public Int16 accsec { get; set; }

        public List<GCElement> elementList { get; set; }

        public GCCell()
        {
            pGCCell();
        }

        void pGCCell()
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

        public void addBox(Int32 x, Int32 y, Int32 b, Int32 h, Int32 layer, Int32 datatype)
        {
            pAddBox(x, y, b, h, layer, datatype);
        }

        void pAddBox(Int32 x, Int32 y, Int32 b, Int32 h, Int32 layer, Int32 datatype)
        {
            addElement(new GCBox(x, y, b, h, layer, datatype));
        }

        public void addBox(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            pAddBox(points, layer, datatype);
        }

        void pAddBox(GeoLibPoint[] points, Int32 layer, Int32 datatype)
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
            Int32 x1, x2, y1, y2;
            p = points[0];
            x1 = p.X;
            y1 = p.Y;
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
                if ((p.X != x1) && (p.X != x2))
                    b = false;
                if ((p.Y != y1) && (p.Y != y2))
                    b = false;
            }
            if (b)
            {
                e = new GCBox(x1, y1, (x2 - x1) + 1, (y2 - y1) + 1, layer, datatype);
                pAddElement(e);
                return;
            }
            else
            {
                e = new GCPolygon(points, layer, datatype);
                pAddElement(e);
                return;
            }
        }

        public GCElement addCircle(Int32 layer, Int32 datatype, GeoLibPoint center, double radius)
        {
            return pAddCircle(layer, datatype, center, radius);
        }

        GCElement pAddCircle(Int32 layer, Int32 datatype, GeoLibPoint center, double radius)
        {
            GCElement e = new GCElement();
            GeoLibPoint[] points = e.ellipse(center, radius, GCSetup.circularDefault);
            e = new GCPolygon(points, layer, datatype);
            return e;
        }

        public void addElement(GCElement element)
        {
            pAddElement(element);
        }

        void pAddElement(GCElement element)
        {
            elementList.Add(element);
        }

        public void addPath(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            pAddPath(points, layer, datatype);
        }

        void pAddPath(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            GCElement e = new GCPath(points, layer, datatype);
            pAddElement(e);
        }

        public void addPolygon(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            pAddPolygon(points, layer, datatype);
        }

        void pAddPolygon(GeoLibPoint[] points, Int32 layer, Int32 datatype)
        {
            GCElement e = new GCPolygon(points, layer, datatype);
            pAddElement(e);
        }

        public void addCellref(GCCell c, GeoLibPoint pos)
        {
            pAddCellref(c, pos);
        }

        void pAddCellref(GCCell c, GeoLibPoint pos)
        {
            GCElement e = new GCCellref(c, pos);
            pAddElement(e);
        }

        public void addCellrefArray(GCCell c, GeoLibPoint[] array, Int32 anzx, Int32 anzy)
        {
            pAddCellrefArray(c, array, anzx, anzy);
        }

        void pAddCellrefArray(GCCell c, GeoLibPoint[] array, Int32 anzx, Int32 anzy)
        {
            GCElement e = new GCCellRefArray(c, array, anzx, anzy);
            pAddElement(e);
        }

        public void addCellrefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, Int32 anzx, Int32 anzy)
        {
            pAddCellrefArray(c, pos1, pos2, anzx, anzy);
        }

        void pAddCellrefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, Int32 anzx, Int32 anzy)
        {
            GCElement e = new GCCellRefArray(c, pos1, pos2, anzx, anzy);
            pAddElement(e);
        }

        public void addCellrefs(int count)
        {
            pAddCellrefs(count);
        }

        void pAddCellrefs(int count)
        {
            GCElement[] e = new GCCellref[count];
            elementList.AddRange(e.ToList());
        }

        public void addCellref()
        {
            pAddCellref();
        }

        void pAddCellref()
        {
            GCElement e = new GCCellref();
            pAddElement(e);
        }

        public void addText(Int32 layer, Int32 datatype, GeoLibPoint pos, string text)
        {
            pAddText(layer, datatype, pos, text);
        }

        void pAddText(Int32 layer, Int32 datatype, GeoLibPoint pos, string text)
        {
            GCElement e = new GCTxt(layer, datatype, pos, text);
            pAddElement(e);
        }

        public bool depend(GCCell cell)
        {
            return pDepend(cell);
        }

        bool pDepend(GCCell cell)
        {
            bool b = false;
            GCCell cellhelp;
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    cellhelp = elementList[f].depend();
                    if (cellhelp == cell)
                    {
                        b = true;
                    }
                    else
                    {
                        if ((cellhelp != null) && (!b))
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

        void pMove(GeoLibPoint p)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].move(p);
                }
            }
        }

        public void moveSelect(GeoLibPoint p)
        {
            pMoveSelect(p);
        }

        void pMoveSelect(GeoLibPoint p)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveSelect(p);
                }
            }
        }

        public void moveToDataType(Int32 datatype)
        {
            pMoveToDataType(datatype);
        }

        void pMoveToDataType(Int32 datatype)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToDataType(datatype);
                }
            }
        }

        public void moveToLayer(Int32 layer)
        {
            pMoveToLayer(layer);
        }

        void pMoveToLayer(Int32 layer)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToLayer(layer);
                }
            }
        }

        public void moveToDataTypeSelect(Int32 datatype)
        {
            pMoveToDataTypeSelect(datatype);
        }

        void pMoveToDataTypeSelect(Int32 datatype)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToDataTypeSelect(datatype);
                }
            }
        }

        public void moveToLayerSelect(Int32 layer)
        {
            pMoveToLayerSelect(layer);
        }

        void pMoveToLayerSelect(Int32 layer)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].moveToLayerSelect(layer);
                }
            }
        }

        public void resize(double size)
        {
            pResize(size);
        }

        void pResize(double size)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].resize(size);
                }
            }
        }

        public void rotateSelect(double angle, GeoLibPoint pos)
        {
            pRotateSelect(angle, pos);
        }

        void pRotateSelect(double angle, GeoLibPoint pos)
        {
            GCStrans m = new GCStrans();
            m.translate(pos.X, pos.Y);
            m.rotate(angle);
            m.translate(-pos.X, -pos.Y);
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if ((elementList[f].isBox()) && (elementList[f].select))
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
        }

        public void scaleSelect(GeoLibPoint pos, GeoLibPoint p2, GeoLibPoint p3)
        {
            pScaleSelect(pos, p2, p3);
        }

        void pScaleSelect(GeoLibPoint pos, GeoLibPoint p2, GeoLibPoint p3)
        {
            GCStrans m = new GCStrans();
            m.translate(pos.X, pos.Y);
            double x = (double)(p3.X - pos.X) / (p2.X - pos.X);
            double y = (double)(p3.Y - pos.Y) / (p2.Y - pos.Y);
            if (x == 0)
            {
                x = 1;
            }
            if (y == 0)
            {
                y = 1;
            }
            if ((p2.X - pos.X) == 0)
            {
                x = 1;
            }
            if ((p2.Y - pos.Y) == 0)
            {
                y = 1;
            }
            m.scale(x, y);
            m.translate(-pos.X, -pos.Y);
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if ((elementList[f].isBox()) && (elementList[f].select))
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
        }

        public void mirrorSelect(GeoLibPoint p1, GeoLibPoint p2)
        {
            pMirrorSelect(p1, p2);
        }

        void pMirrorSelect(GeoLibPoint p1, GeoLibPoint p2)
        {
            GCStrans m = new GCStrans();
            if (p1 == p2)
            {
                m.translate(p1.X, p1.Y);
                m.scale((-1));
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
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if ((elementList[f].isBox()) && (elementList[f].select))
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
        }

        public void flatSelect()
        {
            pFlatSelect();
        }

        void pFlatSelect()
        {
            List<GCElement> e = null;
            List<GCElement> a = null;
            List<GCElement> c = null;
            a = null;
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    e = elementList[f].flatSelect();
                    if (e != null)
                    {
                        c.AddRange(a.ToList());
                        a = e;
                        elementList[f] = null;
                    }
                }
            }
            if (elementList != null)
            {
                elementList.AddRange(a.ToList());
            }
            clean();
        }

        public void cutSelect(GeoLibPoint p1, GeoLibPoint p2)
        {
            pCutSelect(p1, p2);
        }

        void pCutSelect(GeoLibPoint p1, GeoLibPoint p2)
        {
            List<GCElement> a = null;
            List<GCElement> d = null;
            List<GCElement> c = null;
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if ((elementList[f].isBox()) && (elementList[f].select))
                    {
                        GCPolygon p = elementList[f].convertToPolygons()[0];
                        p.select = true;
                        c = p.cutSelect(p1, p2);
                        GCBox b = p.convertToBox();
                        if (b == null)
                        {
                            elementList[f] = p;
                        }
                        else
                        {
                            b.select = true;
                            elementList[f] = b;
                        }
                    }
                    else
                    {
                        c = elementList[f].cutSelect(p1, p2);
                    }
                }
                if (c != null)
                {
                    d.AddRange(a.ToList());
                    a = c;
                }
            }
            elementList.AddRange(a);
            clean();
        }

        public void group(GCCell cell_)
        {
            pGroup(cell_);
        }

        void pGroup(GCCell cell_)
        {
            List<GCElement> a = cell_.elementList.ToList();
            List<GCElement> b;
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].select)
                    {
                        b = new List<GCElement>();
                        b.AddRange(a.ToList());
                        b.Add(elementList[f]);
                        a = b;
                        elementList[f] = null;
                    }
                }
            }
            clean();
            cell_.elementList = a;
        }

        public void minimum(GeoLibPoint pos)
        {
            pMinimum(pos);
        }

        void pMinimum(GeoLibPoint pos)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].minimum(pos);
                }
            }
        }

        public void maximum(GeoLibPoint pos)
        {
            pMaximum(pos);
        }

        void pMaximum(GeoLibPoint pos)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].maximum(pos);
                }
            }
        }

        public void clean()
        {
            pClean();
        }

        void pClean()
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
            elementList.RemoveAll(null);
        }

        public void updateCellref(GCCell oldCell, GCCell newCell)
        {
            pUpdateCellref(oldCell, newCell);
        }

        void pUpdateCellref(GCCell oldCell, GCCell newCell)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    if (elementList[f].depend() == oldCell)
                    {
                        elementList[f].setCellRef(newCell);
                    }
                }
            }
        }

        public bool dependNotSaved()
        {
            return pDependNotSaved();
        }

        bool pDependNotSaved()
        {
            bool b = false;

            if (saved)
            {
                return false;
            }

            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    GCCell cellHelp = elementList[f].depend();
                    if ((cellHelp != null) && (!b))
                    {
                        if (!cellHelp.saved)
                        {
                            b = true;
                        }
                    }
                }
            }

            return b;
        }

        public void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            //bgnstr
            gw.bw.Write((UInt16)28);
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
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].saveGDS(gw);
                }
            }
            //endstr
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)7);
            gw.bw.Write((byte)0);
            saved = true;
        }

        public void saveOASIS(oasWriter ow)
        {
            pSaveOASIS(ow);
        }

        void pSaveOASIS(oasWriter ow)
        {
            for (int f = 0; f < elementList.Count; f++)
            {
                if (elementList[f] != null)
                {
                    elementList[f].saveOASIS(ow);
                }
            }
            saved = true;
        }
    }
}
