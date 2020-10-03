using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;

namespace geoCoreLib
{
    public class GCBox : GCElement
    {
        GeoLibRectangle rect;

        public GCBox()
        {
        }

        public GCBox(GeoLibRectangle r, Int32 layer, Int32 datatype)
        {
            pGCBox(r, layer, datatype);
        }

        void pGCBox(GeoLibRectangle r, Int32 layer, Int32 datatype)
        {
            rect = r;
            layer_nr = layer;
            datatype_nr = datatype;
        }

        public GCBox(Int32 x_, Int32 y_, Int32 b_, Int32 h_, Int32 layer, Int32 datatype)
        {
            pGCBox(x_, y_, b_, h_, layer, datatype);
        }

        void pGCBox(Int32 x_, Int32 y_, Int32 b_, Int32 h_, Int32 layer, Int32 datatype)
        {
            rect = new GeoLibRectangle(x_, y_, b_, h_);
            pGCBox(rect, layer, datatype);
        }

        public override void clean()
        {
            // Nothing needed since the rectangle primitive does everything.
        }

        public override bool correct()
        {
            return pCorrect();
        }

        bool pCorrect()
        {
            if (rect.Left == rect.Right)
            {
                return false;
            }
            if (rect.Top == rect.Bottom)
            {
                return false;
            }
            return true;
        }

        public override void maximum(GeoLibPoint p)
        {
            pMaximum(p);
        }

        void pMaximum(GeoLibPoint p)
        {
            p.X = rect.Right;
            p.Y = rect.Top;
        }

        public override void minimum(GeoLibPoint p)
        {
            pMinimum(p);
        }

        void pMinimum(GeoLibPoint p)
        {
            p.X = rect.Left;
            p.Y = rect.Bottom;
        }

        public override void move(GeoLibPoint pos)
        {
            pMove(pos);
        }

        void pMove(GeoLibPoint pos)
        {
            rect.Offset(pos);
        }

        public override void moveSelect(GeoLibPoint p)
        {
            pMoveSelect(p);
        }

        void pMoveSelect(GeoLibPoint p)
        {
            if (select)
            {
                rect.Offset(p);
            }
        }

        public override void resize(double size)
        {
            pResize(size);
        }

        void pResize(double size)
        {
            rect = new GeoLibRectangle(
                    rect.X,
                    rect.Y,
                    (Int32)(rect.Width * size),
                    (Int32)(rect.Height * size)
                );
        }

        public override List<GCPolygon> convertToPolygons()
        {
            return pConvertToPolygons();
        }

        List<GCPolygon> pConvertToPolygons()
        {
            List<GCPolygon> ret = new List<GCPolygon>();
            GeoLibPoint[] points = new GeoLibPoint[5];
            points[0] = new GeoLibPoint(rect.Left, rect.Top);
            points[1] = new GeoLibPoint(rect.Right, rect.Top);
            points[2] = new GeoLibPoint(rect.Right, rect.Bottom);
            points[3] = new GeoLibPoint(rect.Left, rect.Bottom);
            points[4] = new GeoLibPoint(rect.Left, rect.Top);
            ret.Add(new GCPolygon(points, layer_nr, datatype_nr));
            return ret;
        }

        public override void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            // box
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x2D);
            gw.bw.Write((byte)0);
            //layer
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x0D);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)layer_nr);
            //boxtype
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x2E);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)datatype_nr);
            //xy 
            gw.bw.Write((UInt16)((5 * 2 * 4) + 4));
            gw.bw.Write((byte)0x10);
            gw.bw.Write((byte)3);
            gw.bw.Write(rect.Left);
            gw.bw.Write(rect.Top);
            gw.bw.Write(rect.Right);
            gw.bw.Write(rect.Top);
            gw.bw.Write(rect.Right);
            gw.bw.Write(rect.Bottom);
            gw.bw.Write(rect.Left);
            gw.bw.Write(rect.Bottom);
            gw.bw.Write(rect.Left);
            gw.bw.Write(rect.Top);
            // endel
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x11);
            gw.bw.Write((byte)0);
        }

        public override void saveOASIS(oasWriter ow)
        {
            pSaveOASIS(ow);
        }

        void pSaveOASIS(oasWriter ow)
        {
            if (!ow.modal.absoluteMode)
            {
                ow.modal.absoluteMode = true;
            }

            byte info_byte = 0;
            if (layer_nr != ow.modal.layer)
            {
                info_byte += 1;
            }
            if (datatype_nr != ow.modal.datatype)
            {
                info_byte += 2;
            }
            if (rect.Left != ow.modal.geometry_x)
            {
                info_byte += 16;
            }
            if (rect.Top != ow.modal.geometry_y)
            {
                info_byte += 8;
            }
            if (rect.Bottom - rect.Top != ow.modal.geometry_h)
            {
                info_byte += 32;
            }
            if (rect.Right - rect.Left != ow.modal.geometry_w)
            {
                info_byte += 64;
            }
            if (rect.Right - rect.Left == rect.Bottom - rect.Top)
            {
                info_byte += 128;
                info_byte = (byte)(info_byte & (255 - 32));
            }
            ow.writeUnsignedInteger(20);
            ow.writeRaw(info_byte);
            if ((info_byte & 1) > 0)
            {
                ow.modal.layer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.layer);
            }
            if ((info_byte & 2) > 0)
            {
                ow.modal.datatype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
            }
            if ((info_byte & 64) > 0)
            {
                ow.modal.geometry_w = rect.Right - rect.Left;
                ow.writeUnsignedInteger((uint)ow.modal.geometry_w);
                if ((info_byte & 128) > 0)
                {
                    ow.modal.geometry_h = ow.modal.geometry_w;
                }
            }
            if ((info_byte & 32) > 0)
            {
                ow.modal.geometry_h = rect.Bottom - rect.Top;
                ow.writeUnsignedInteger((uint)ow.modal.geometry_h);
            }
            if ((info_byte & 16) > 0)
            {
                ow.modal.geometry_x = rect.Left;
                ow.writeSignedInteger(ow.modal.geometry_x);
            }
            if ((info_byte & 8) > 0)
            {
                ow.modal.geometry_y = rect.Top;
                ow.writeSignedInteger(ow.modal.geometry_y);
            }
        }
    }
}
