using gds;
using geoLib;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoCoreLib
{
    public class GCPath : GCElement
    {
        GeoLibPoint[] pointarray;
        int width;
        int cap;
        bool isRound;

        public GCPath()
        {
            pGCPath();
        }

        void pGCPath()
        {
            pointarray = new GeoLibPoint[0];
        }

        public GCPath(GeoLibPoint[] points, int layer, int datatype)
        {
            pGCPath(points, layer, datatype);
        }

        void pGCPath(GeoLibPoint[] points, int layer, int datatype)
        {
            layer_nr = layer;
            datatype_nr = datatype;
            pointarray = points;
            width = 0;
            cap = 0;
            clean();
        }

        public override bool correct()
        {
            return pCorrect();
        }

        bool pCorrect()
        {
            if (pointarray.Count() < 2)
            {
                return false;
            }
            if ((pointarray.Count() == 2) && (pointarray[0] == pointarray[1]) && (cap == 0))
            {
                return false;
            }
            return true;
        }

        public override void minimum(GeoLibPoint pos)
        {
            pMinimum(pos);
        }

        void pMinimum(GeoLibPoint pos)
        {
            GeoLibPoint p = GeoWrangler.getMinimumPoint(pointarray);
            pos.X = Math.Min(pos.X, p.X);
            pos.Y = Math.Min(pos.Y, p.Y);
        }

        public override void maximum(GeoLibPoint pos)
        {
            pMaximum(pos);
        }

        void pMaximum(GeoLibPoint pos)
        {
            GeoLibPoint p = GeoWrangler.getMaximumPoint(pointarray);
            pos.X = Math.Max(pos.X, p.X);
            pos.Y = Math.Max(pos.Y, p.Y);
        }

        public override void moveSelect(GeoLibPoint pos)
        {
            pMoveSelect(pos);
        }

        void pMoveSelect(GeoLibPoint pos)
        {
            if (select)
            {
                int pointArrayCount = pointarray.Length;
#if GCTHREADED
                Parallel.For(i, pointArrayCount, (i) =>
#else
                for (int i = 0; i < pointArrayCount; i++)
#endif
                {
                    pointarray[i].Offset(pos);
                }
            }
#if GCTHREADED
            );
#endif
        }

        public override void move(GeoLibPoint pos)
        {
            pMove(pos);
        }

        void pMove(GeoLibPoint pos)
        {
            int pointArrayCount = pointarray.Length;
#if GCTHREADED
            Parallel.For(i, pointArrayCount, (i) =>
#else
            for (int i = 0; i < pointArrayCount; i++)
#endif
            {
                pointarray[i].Offset(pos);
            }
#if GCTHREADED
            );
#endif
        }

        public override void resize(double size)
        {
            pResize(size);
        }

        void pResize(double size)
        {
            width = (Int32)(size * width);
            pointarray = GeoWrangler.resize(pointarray, size);
        }

        public override void clean()
        {
            pClean();
        }

        void pClean()
        {
            for (int i = 0; i < pointarray.Length - 1; i++)
            {
                if (pointarray.Length <= 2)
                    return; //no area
                if (pointarray[i] == pointarray[i + 1])
                {
                    deletePoint(i + 1);
                    i--;
                }
            }
        }

        public void deletePoint(int pos)
        {
            pDeletePoint(pos);
        }

        void pDeletePoint(int pos)
        {
            List<GeoLibPoint> newArray = pointarray.ToList();
            newArray.RemoveAt(pos);
            pointarray = newArray.ToArray();
        }

        public override void setWidth(int width)
        {
            pSetWidth(width);
        }

        void pSetWidth(int width)
        {
            this.width = width;
        }

        public override void setCap(int cap)
        {
            pSetCap(cap);
        }

        void pSetCap(int cap)
        {
            this.cap = cap;
        }

        public override List<GCPolygon> convertToPolygons()
        {
            return pConvertToPolygons();
        }

        List<GCPolygon> pConvertToPolygons()
        {
            List<GCPolygon> ret = new List<GCPolygon>();
            GeoLibPoint[] tmp = GeoWrangler.inflatePath(pointarray, width);
            ret.Add(new GCPolygon(tmp, layer_nr, datatype_nr));
            return ret;
        }

        public void setRound(bool val)
        {
            pSetRound(val);
        }

        void pSetRound(bool val)
        {
            isRound = val;
        }

        public override void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            // path 
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)9);
            gw.bw.Write((byte)0);
            //layer
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x0D);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)layer_nr);
            //datatype
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x0E);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)(datatype_nr));
            //pathtype
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x21);
            gw.bw.Write((byte)2);
            if (!isRound)
            {
                gw.bw.Write((Int16)cap);
            }
            else
            {
                gw.bw.Write((Int16)1);
            }
            //width
            gw.bw.Write((UInt16)8);
            gw.bw.Write((byte)0x0F);
            gw.bw.Write((byte)3);
            gw.bw.Write(width);

            int i = pointarray.Length;
            if (i > 8191)
            {
                i = 8191;
                // "Path with more than 8191 points. Data is lost."
            }
            //xy 
            int val = (i * 2 * 4) + 4;
            gw.bw.Write((UInt16)val);
            gw.bw.Write((byte)0x10);
            gw.bw.Write((byte)3);
            for (int k = 0; k < i; k++)
            {
                gw.bw.Write(pointarray[k].X);
                gw.bw.Write(pointarray[k].Y);
            }
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
            byte info_byte = 128 + 32;  //write pointlist+extension;

            if (layer_nr != ow.modal.layer)
            {
                info_byte += 1;
            }

            if (datatype_nr != ow.modal.datatype)
            {
                info_byte += 2;
            }

            if (pointarray[0].X != ow.modal.geometry_x)
            {
                info_byte += 16;
            }

            if (pointarray[0].Y != ow.modal.geometry_y)
            {
                info_byte += 8;
            }

            // Yes, indeed. Oasis is limited to path widths that are factors of 2. Strange, but there you are.
            if ((width / 2) != ow.modal.geometry_w)
            {
                info_byte += 64;
            }

            ow.writeUnsignedInteger(22);
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
                // Yes, indeed. Oasis is limited to path widths that are factors of 2. Strange, but there you are.
                ow.modal.geometry_w = width / 2;
                if (ow.modal.geometry_w < 0)
                {
                    ow.modal.geometry_w = -ow.modal.geometry_w;
                }
                ow.writeUnsignedInteger((uint)ow.modal.geometry_w);
            }
            // 2 means the end gets extended by the half the path width.
            if (cap == 0)
            {
                ow.writeUnsignedInteger(5);
            }
            else if (cap == 2)
            {
                ow.writeUnsignedInteger(10);
            }
            else
            {
                ow.writeUnsignedInteger(5);
                // Round caps are not possible in OASIS files.;
            }
            ow.writePointArray(pointarray, false);
            if ((info_byte & 16) > 0)
            {
                ow.modal.geometry_x = pointarray[0].X;
                ow.writeSignedInteger(ow.modal.geometry_x);
            }
            if ((info_byte & 8) > 0)
            {
                ow.modal.geometry_y = pointarray[0].Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
            }
        }
    }
}
