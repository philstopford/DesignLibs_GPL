using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using utility;

namespace geoCoreLib
{
    public class GCCellref : GCElement
    {
        public GeoLibPoint point { get; set; }
        string name;
        public GCCell cell_ref { get; set; }
        public GCStrans trans { get; set; }

        public override void setName(string s)
        {
            pSetName(s);
        }

        void pSetName(string s)
        {
            name = s;
        }

        public GCCellref(GCCell c, GeoLibPoint pos)
        {
            pGCCellref(c, pos);
        }

        void pGCCellref(GCCell c, GeoLibPoint pos)
        {
            cell_ref = c;
            // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
            layer_nr = -1;
            datatype_nr = -1;
            point = pos;
            trans = new GCStrans();
            trans.reset();
        }

        public GCCellref()
        {
            pGCCellref();
        }

        void pGCCellref()
        {
            cell_ref = null;
            // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
            layer_nr = -1;
            datatype_nr = -1;
            point = new GeoLibPoint(0, 0);
            trans = new GCStrans();
            trans.reset();
        }

        public override void setPos(GeoLibPoint p)
        {
            pSetPos(p);
        }

        void pSetPos(GeoLibPoint p)
        {
            point = new GeoLibPoint(p.X, p.Y);
        }

        public override GeoLibPoint getPos()
        {
            return pGetPos();
        }

        GeoLibPoint pGetPos()
        {
            return point;
        }

        public override void minimum(GeoLibPoint p)
        {
            pMinimum(p);
        }

        void pMinimum(GeoLibPoint p)
        {
            GeoLibPoint pos1, pos2;
            pos1 = new GeoLibPoint(p.X - point.X, p.Y - point.Y);
            if (trans.mirror_x)
            {
                pos1.Y = -pos1.Y;
            }
            if (trans.angle == 180)
            {
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
            }
            pos2 = pos1;
            cell_ref.maximum(pos1);
            cell_ref.minimum(pos2);
            if (trans.mirror_x)
            {
                pos1.Y = -pos1.Y;
                pos2.Y = -pos2.Y;
            }
            if (trans.angle == 180)
            {
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
            }
            if (trans.angle == 180)
            {
                pos2.Y = -pos2.Y;
                pos2.X = -pos2.X;
            }
            pos1.Offset(point);
            pos2.Offset(point);
            if (pos2.X < p.X)
            {
                p.X = pos2.X;
            }
            if (pos1.X < p.X)
            {
                p.X = pos1.X;
            }
            if (pos2.Y < p.Y)
            {
                p.Y = pos2.Y;
            }
            if (pos1.Y < p.Y)
            {
                p.Y = pos1.Y;
            }
        }

        public override void scale(double factor)
        {
            pScale(factor);
        }

        void pScale(double factor)
        {
            trans.mag = factor;
        }

        public override void rotate(double angle)
        {
            pRotate(angle);
        }

        void pRotate(double angle)
        {
            trans.angle = angle;
        }

        public override void maximum(GeoLibPoint p)
        {
            pMaximum(p);
        }

        void pMaximum(GeoLibPoint p)
        {
            GeoLibPoint pos1, pos2;
            pos1 = new GeoLibPoint(p.X - point.X, p.Y - point.Y);
            if (trans.mirror_x)
            {
                pos1.Y = -pos1.Y;
            }
            if (trans.angle == 180)
            {
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
            }
            pos2 = pos1;
            cell_ref.maximum(pos1);
            cell_ref.minimum(pos2);
            if (trans.mirror_x)
            {
                pos1.Y = -pos1.Y;
                pos2.Y = -pos2.Y;
            }
            if (trans.angle == 180)
            {
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
            }
            if (trans.angle == 180)
            {
                pos2.Y = -pos2.Y;
                pos2.X = -pos2.X;
            }
            pos1.Offset(point);
            pos2.Offset(point);
            if (pos2.X > p.X)
            {
                p.X = pos2.X;
            }
            if (pos1.X > p.X)
            {
                p.X = pos1.X;
            }
            if (pos2.Y > p.Y)
            {
                p.Y = pos2.Y;
            }
            if (pos1.Y > p.Y)
            {
                p.Y = pos1.Y;
            }
        }

        public override GCCell depend()
        {
            return pDepend();
        }

        GCCell pDepend()
        {
            return cell_ref;
        }

        public override void move(GeoLibPoint p)
        {
            pMove(p);
        }

        void pMove(GeoLibPoint p)
        {
            point.Offset(p);
        }

        public override void moveSelect(GeoLibPoint p)
        {
            pMoveSelect(p);
        }

        void pMoveSelect(GeoLibPoint p)
        {
            point.Offset(p);
        }

        public override void setMirrorx()
        {
            pSetMirrors();
        }

        void pSetMirrors()
        {
            trans.setMirror_x();
        }

        public override void resize(double factor)
        {
            pResize(factor);
        }

        void pResize(double factor)
        {
            point.X = (Int32)(point.X * factor);
            point.Y = (Int32)(point.Y * factor);
        }

        public override bool isCellref()
        {
            return pIsCellref();
        }

        bool pIsCellref()
        {
            return true;
        }

        public override void setCellRef(GCCell cellRef)
        {
            pSetCellRef(cellRef);
        }

        void pSetCellRef(GCCell cellRef)
        {
            cell_ref = cellRef;
        }

        public override GCCell getCellref()
        {
            return pGetCellRef();
        }

        GCCell pGetCellRef()
        {
            return cell_ref;
        }

        public override void saveGDS(gdsWriter bw)
        {
            pSaveGDS(bw);
        }

        void pSaveGDS(gdsWriter gw)
        {
            // Guard against incoming broken cellref information
            if (cell_ref == null)
            {
                return;
            }
            //SRef
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x0A);
            gw.bw.Write((byte)0);
            gw.writeString(cell_ref.cellName, 0x12);
            int strans_ = 0;
            if (trans.mirror_x)
            {
                strans_ |= 32768;
            }
            //STRANS
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x1A);
            gw.bw.Write((byte)1);
            gw.bw.Write((Int16)strans_);
            //mag
            gw.bw.Write((UInt16)12);
            gw.bw.Write((byte)0x1B);
            gw.bw.Write((byte)5);
            gw.write8ByteReal(trans.mag);
            //angle
            gw.bw.Write((UInt16)12);
            gw.bw.Write((byte)0x1C);
            gw.bw.Write((byte)5);
            if ((trans.mirror_x) && (trans.angle != 0))
            {
                gw.write8ByteReal(360 - trans.angle);
            }
            else
            {
                gw.write8ByteReal(trans.angle);
            }
            //xy
            gw.bw.Write((UInt16)((2 * 4) + 4));
            gw.bw.Write((byte)0x10);
            gw.bw.Write((byte)3);
            gw.bw.Write(point.X);
            gw.bw.Write(point.Y);
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
            byte info_byte = 128;  //explicid cellname;
            if (trans.mirror_x)
            {
                info_byte += 1;
            }
            if (trans.mag != 1)
            {
                info_byte += 4;
            }
            if (trans.angle != 0)
            {
                info_byte += 2;
            }
            if (point.X != ow.modal.placement_x)
            {
                info_byte += 32;
            }
            if (point.Y != ow.modal.placement_y)
            {
                info_byte += 16;
            }
            ow.writeUnsignedInteger(18);
            ow.writeRaw(info_byte);
            ow.writeString(cell_ref.cellName);
            if ((info_byte & 4) > 0)
            {
                ow.writeReal(trans.mag);
            }
            if ((info_byte & 2) > 0)
            {
                ow.writeReal(trans.angle);
            }
            if ((info_byte & 32) > 0)
            {
                ow.modal.placement_x = point.X;
                ow.writeSignedInteger(ow.modal.placement_x);
            }
            if ((info_byte & 16) > 0)
            {
                ow.modal.placement_y = point.Y;
                ow.writeSignedInteger(ow.modal.placement_y);
            }
        }

        public override List<GCPolygon> convertToPolygons()
        {
            return pConvertToPolygons();
        }

        List<GCPolygon> pConvertToPolygons()
        {
            List<GCPolygon> ret = new List<GCPolygon>();

            for (int element = 0; element < cell_ref.elementList.Count; element++)
            {
                List<GCPolygon> tmp = cell_ref.elementList[element].convertToPolygons();
                ret.AddRange(tmp);
            }

#if GCTHREADED
            Parallel.For(0, ret.Count, (poly, loopstate) =>
#else
            for (int poly = 0; poly < ret.Count; poly++)
#endif
            {
                ret[poly].rotate(trans.angle, point);
                ret[poly].scale(trans.mag);
            }
#if GCTHREADED
            );
#endif
            return ret;
        }
    }
}
