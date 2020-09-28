using gds;
using geoLib;
using oasis;
using System;

namespace geoCoreLib
{
    public class GCCellRefArray : GCElement
    {
        GeoLibPoint point;
        public Int32 count_x, count_y;
        public GeoLibPoint space;
        public GCCell cell_ref { get; set; }
        string name;
        public GCStrans trans { get; set; }

        public GCCellRefArray(GCCell c, GeoLibPoint[] array, Int32 xCount, Int32 yCount)
        {
            cell_ref = c;
            point = array[0];
            space = new GeoLibPoint((array[1].X - point.X) / xCount,  (array[2].Y - point.Y) / yCount);
            count_x = xCount;
            count_y = yCount;
            // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
            layer_nr = -1;
            datatype_nr = -1;
            trans = new GCStrans();
            trans.reset();
        }

        public GCCellRefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, Int32 xCount, Int32 yCount)
        {
            cell_ref = c;
            point = pos1;
            GeoLibPoint p = new GeoLibPoint(pos2.X - pos1.X, pos2.Y - pos1.Y);
            space = new GeoLibPoint(p.X, p.Y);
            count_x = xCount;
            count_y = yCount;
            // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
            layer_nr = -1;
            datatype_nr = -1;
            trans = new GCStrans();
            trans.reset();
        }

        public GCCellRefArray()
        {
            cell_ref = null;
            // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
            layer_nr = -1;
            datatype_nr = -1;
            trans = new GCStrans();
            trans.reset();
        }

        public override GCCell depend()
        {
            return pDepend();
        }

        GCCell pDepend()
        {
            return cell_ref;
        }

        public override void minimum(GeoLibPoint p)
        {
            pMinimum(p);
        }

        void pMinimum(GeoLibPoint p)
        {
            GeoLibPoint pos1, pos2, pos3;
            for (int x = 0; x < 2; x++)
            {
                for (int y = 0; y < 2; y++)
                {
                    pos3 = point;
                    pos3.Offset(new GeoLibPoint(space.X * x * (count_x - 1), space.Y * y * (count_y - 1)));
                    pos1 = new GeoLibPoint(p.X - pos3.X, p.Y - pos3.Y);
                    if (trans.mirror_x)
                    {
                        pos1.Y = (-pos1.Y);
                    };
                    pos2 = pos1;
                    cell_ref.maximum(pos1);
                    cell_ref.minimum(pos2);
                    if (trans.mirror_x)
                    {
                        pos1.Y = (-pos1.Y);
                        pos2.Y = (-pos2.Y);
                    };
                    pos1.Offset(pos3);
                    pos2.Offset(pos3);
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
            }
        }

        public override void maximum(GeoLibPoint p)
        {
            pMaximum(p);
        }

        void pMaximum(GeoLibPoint p)
        {
            GeoLibPoint pos1, pos2, pos3;
            for (int x = 0; x < 2; x++)
            {
                for (int y = 0; y < 2; y++)
                {
                    pos3 = point;
                    pos3.Offset(new GeoLibPoint(space.X * x * (count_x - 1), space.Y * y * (count_y - 1)));
                    pos1 = new GeoLibPoint(p.X - pos3.X, p.Y - pos3.Y);
                    if (trans.mirror_x)
                    {
                        pos1.Y = -pos1.Y;
                    };
                    pos2 = pos1;
                    cell_ref.maximum(pos1);
                    cell_ref.minimum(pos2);
                    if (trans.mirror_x)
                    {
                        pos1.Y = -pos1.Y;
                        pos2.Y = -pos2.Y;
                    };
                    pos1.Offset(pos3);
                    pos2.Offset(pos3);
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
            }
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
            if (select)
            {
                point.Offset(p);
            }
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
            space.X = (Int32)(space.X * factor);
            space.Y = (Int32)(space.Y * factor);
        }

        public override void setPos(GeoLibPoint p)
        {
            pSetPos(p);
        }

        void pSetPos(GeoLibPoint p)
        {
            point = p;
        }

        public override void setCellRef(GCCell cellRef)
        {
            pSetCellRef(cellRef);
        }

        void pSetCellRef(GCCell cellRef)
        {
            cell_ref = cellRef;
        }

        public override void saveGDS(gdsWriter gw)
        {
            pSaveGDS(gw);
        }

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        void pSaveGDS(gdsWriter gw)
        {
            //SRef
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x0B);
            gw.bw.Write((byte)0);
            gw.writeString(cell_ref.cellName, 0x12);
            Int32 strans_ = 0;
            if (trans.mirror_x)
            {
                strans_ |= 0x8000;
            };
            //STRANS
            gw.bw.Write((UInt16)6);
            gw.bw.Write((byte)0x1A);
            gw.bw.Write((byte)1);
            gw.bw.Write((Int16)(strans_));
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
            //colrow
            gw.bw.Write((UInt16)8);
            gw.bw.Write((byte)0x13);
            gw.bw.Write((byte)2);
            gw.bw.Write((Int16)count_x);
            gw.bw.Write((Int16)count_y);
            //xy
            gw.bw.Write((UInt16)((3 * 2 * 4) + 4));
            gw.bw.Write((byte)0x10);
            gw.bw.Write((byte)3);
            gw.bw.Write(point.X);
            gw.bw.Write(point.Y);
            GeoLibPoint pos = new GeoLibPoint(space.X * count_x + point.X, space.Y + point.Y);
            gw.bw.Write(pos.X);
            gw.bw.Write(pos.Y);
            pos = new GeoLibPoint(space.X * +point.X, space.Y * count_y + point.Y);
            gw.bw.Write(pos.X);
            gw.bw.Write(pos.Y);
            // endel
            gw.bw.Write((UInt16)4);
            gw.bw.Write((byte)0x11);
            gw.bw.Write((byte)0);
        }

        public override void saveOASIS(oasWriter ow)
        {
            pSaveOASIS(ow);
        }

        // Cellrefarrays also have to resolve to integer placement.
        // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
        void pSaveOASIS(oasWriter ow)
        {
            if (!ow.modal.absoluteMode)
            {
                ow.modal.absoluteMode = true;
            }
            byte info_byte = 128;  //explicid cellname,repition;
            if ((count_x > 1) || (count_y > 1))
            {
                info_byte += 8;
            }

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
            if ((info_byte & 8) > 0)
            {
                if (count_x == 1)
                {
                    ow.writeUnsignedInteger(3);
                    ow.modal.y_dimension = count_y;
                    ow.writeUnsignedInteger((uint)(count_y - 2));
                    ow.modal.y_space = space.Y;
                    ow.writeUnsignedInteger((uint)space.Y);
                }
                else if (count_y == 1)
                {
                    ow.writeUnsignedInteger(2);
                    ow.modal.x_dimension = count_x;
                    ow.writeUnsignedInteger((uint)(count_x - 2));
                    ow.modal.x_space = space.X;
                    ow.writeUnsignedInteger((uint)space.X);
                }
                else
                {
                    ow.writeUnsignedInteger(1);
                    ow.modal.x_dimension = count_x;
                    ow.modal.y_dimension = count_y;
                    ow.writeUnsignedInteger((uint)(count_x - 2));
                    ow.writeUnsignedInteger((uint)(count_y - 2));
                    ow.modal.x_space = space.X;
                    ow.modal.y_space = space.Y;
                    ow.writeUnsignedInteger((uint)space.X);
                    ow.writeUnsignedInteger((uint)space.Y);
                }
            }
        }

        public override bool isCellrefArray()
        {
            return pIsCellRefArray();
        }

        bool pIsCellRefArray()
        {
            return true;
        }

        public override GCCell getCellref() {
            return pGetCellRef();
        }

        GCCell pGetCellRef()
        {
            return cell_ref;
        }

        public override GeoLibPoint getPos()
        {
            return pGetPos();
        }

        GeoLibPoint pGetPos()
        {
            return point;
        }
    }
}
