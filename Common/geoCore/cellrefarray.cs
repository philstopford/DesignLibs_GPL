using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoCoreLib;

public class GCCellRefArray : GCElement
{
    private GeoLibPoint point;
    public int count_x, count_y;
    public GeoLibPoint pitch;
    public GCCell cell_ref { get; set; }
    private string name;
    public GCStrans trans { get; set; }

    public GCCellRefArray(GCCell c, GeoLibPoint[] array, int xCount, int yCount)
    {
        cell_ref = c;
        point = array[0];
        pitch = new GeoLibPoint((array[1].X - point.X) / xCount,  (array[2].Y - point.Y) / yCount);
        count_x = xCount;
        count_y = yCount;
        // Tag layer and datatype to allow this element to be filtered out from LD and geo lists.
        layer_nr = -1;
        datatype_nr = -1;
        trans = new GCStrans();
        trans.reset();
    }

    public GCCellRefArray(GCCell c, GeoLibPoint pos1, GeoLibPoint pos2, int xCount, int yCount)
    {
        cell_ref = c;
        point = pos1;
        GeoLibPoint p = new(pos2.X - pos1.X, pos2.Y - pos1.Y);
        pitch = new GeoLibPoint(p.X, p.Y);
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

    private GCCell pDepend()
    {
        return cell_ref;
    }

    public override void minimum(GeoLibPoint p)
    {
        pMinimum(p);
    }

    private void pMinimum(GeoLibPoint p)
    {
        for (int x = 0; x < 2; x++)
        {
            for (int y = 0; y < 2; y++)
            {
                GeoLibPoint pos3 = point;
                pos3.Offset(new GeoLibPoint(pitch.X * x * (count_x - 1), pitch.Y * y * (count_y - 1)));
                GeoLibPoint pos1 = new(p.X - pos3.X, p.Y - pos3.Y);
                pos1.Y = trans.mirror_x switch
                {
                    true => -pos1.Y,
                    _ => pos1.Y
                };
                GeoLibPoint pos2 = pos1;
                cell_ref.maximum(pos1);
                cell_ref.minimum(pos2);
                switch (trans.mirror_x)
                {
                    case true:
                        pos1.Y = -pos1.Y;
                        pos2.Y = -pos2.Y;
                        break;
                }
                pos1.Offset(pos3);
                pos2.Offset(pos3);
                p.X = Math.Min(p.X, pos2.X);
                p.X = Math.Min(p.X, pos1.X);
                p.Y = Math.Min(p.Y, pos2.Y);
                p.Y = Math.Min(p.Y, pos1.Y);
            }
        }
    }

    public override void maximum(GeoLibPoint p)
    {
        pMaximum(p);
    }

    private void pMaximum(GeoLibPoint p)
    {
        for (int x = 0; x < 2; x++)
        {
            for (int y = 0; y < 2; y++)
            {
                GeoLibPoint pos3 = point;
                pos3.Offset(new GeoLibPoint(pitch.X * x * (count_x - 1), pitch.Y * y * (count_y - 1)));
                GeoLibPoint pos1 = new(p.X - pos3.X, p.Y - pos3.Y);
                pos1.Y = trans.mirror_x switch
                {
                    true => -pos1.Y,
                    _ => pos1.Y
                };
                GeoLibPoint pos2 = pos1;
                cell_ref.maximum(pos1);
                cell_ref.minimum(pos2);
                switch (trans.mirror_x)
                {
                    case true:
                        pos1.Y = -pos1.Y;
                        pos2.Y = -pos2.Y;
                        break;
                }
                pos1.Offset(pos3);
                pos2.Offset(pos3);
                p.X = Math.Max(p.X, pos1.X);
                p.X = Math.Max(p.X, pos2.X);
                p.Y = Math.Max(p.Y, pos1.Y);
                p.Y = Math.Max(p.Y, pos2.Y);
            }
        }
    }

    public override void move(GeoLibPoint p)
    {
        pMove(p);
    }

    private void pMove(GeoLibPoint p)
    {
        point.Offset(p);
    }

    public override void moveSelect(GeoLibPoint p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(GeoLibPoint p)
    {
        switch (select)
        {
            case true:
                point.Offset(p);
                break;
        }
    }

    public override void setMirrorx()
    {
        pSetMirrors();
    }

    private void pSetMirrors()
    {
        trans.setMirror_x();
    }

    public override void resize(double factor)
    {
        pResize(factor);
    }

    private void pResize(double factor)
    {
        point.X = (int)(point.X * factor);
        point.Y = (int)(point.Y * factor);
        pitch.X = (int)(pitch.X * factor);
        pitch.Y = (int)(pitch.Y * factor);
    }

    public override void setPos(GeoLibPoint p)
    {
        pSetPos(p);
    }

    private void pSetPos(GeoLibPoint p)
    {
        point = p;
    }

    public override void setCellRef(GCCell cellRef)
    {
        pSetCellRef(cellRef);
    }

    private void pSetCellRef(GCCell cellRef)
    {
        cell_ref = cellRef;
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    // Cellrefarrays also have to resolve to integer placement.
    // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
    private void pSaveGDS(gdsWriter gw)
    {
        //SRef
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x0B);
        gw.bw.Write((byte)0);
        gw.writeString(cell_ref.cellName, 0x12);
        int strans_ = 0;
        switch (trans.mirror_x)
        {
            case true:
                strans_ |= 0x8000;
                break;
        }
        //STRANS
        gw.bw.Write((ushort)6);
        gw.bw.Write((byte)0x1A);
        gw.bw.Write((byte)1);
        gw.bw.Write((short)strans_);
        //mag
        gw.bw.Write((ushort)12);
        gw.bw.Write((byte)0x1B);
        gw.bw.Write((byte)5);
        gw.write8ByteReal(trans.mag);
        //angle
        gw.bw.Write((ushort)12);
        gw.bw.Write((byte)0x1C);
        gw.bw.Write((byte)5);
        switch (trans.mirror_x)
        {
            case true when trans.angle != 0:
                gw.write8ByteReal(360 - trans.angle);
                break;
            default:
                gw.write8ByteReal(trans.angle);
                break;
        }
        //colrow
        gw.bw.Write((ushort)8);
        gw.bw.Write((byte)0x13);
        gw.bw.Write((byte)2);
        gw.bw.Write((short)count_x);
        gw.bw.Write((short)count_y);
        //xy
        gw.bw.Write((ushort)(3 * 2 * 4 + 4));
        gw.bw.Write((byte)0x10);
        gw.bw.Write((byte)3);
        gw.bw.Write(point.X);
        gw.bw.Write(point.Y);
        GeoLibPoint pos = new(pitch.X * count_x + point.X, pitch.Y + point.Y);
        gw.bw.Write(pos.X);
        gw.bw.Write(pos.Y);
        pos = new GeoLibPoint(pitch.X * +point.X, pitch.Y * count_y + point.Y);
        gw.bw.Write(pos.X);
        gw.bw.Write(pos.Y);
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x11);
        gw.bw.Write((byte)0);
    }

    public override void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    // Cellrefarrays also have to resolve to integer placement.
    // Placement errors will occur if the x, y instance counts do not divide the array X, Y values cleanly.
    private void pSaveOASIS(oasWriter ow)
    {
        ow.modal.absoluteMode = ow.modal.absoluteMode switch
        {
            false => true,
            _ => ow.modal.absoluteMode
        };
        byte info_byte = 128;  //explicid cellname,repition;
        if (count_x > 1 || count_y > 1)
        {
            info_byte += 8;
        }

        switch (trans.mirror_x)
        {
            case true:
                info_byte += 1;
                break;
        }

        switch (Math.Abs(trans.mag - 1))
        {
            case > double.Epsilon:
                info_byte += 4;
                break;
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

        switch (info_byte & 4)
        {
            case > 0:
                ow.writeReal(trans.mag);
                break;
        }

        switch (info_byte & 2)
        {
            case > 0:
                ow.writeReal(trans.angle);
                break;
        }
        switch (info_byte & 32)
        {
            case > 0:
                ow.modal.placement_x = point.X;
                ow.writeSignedInteger(ow.modal.placement_x);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.placement_y = point.Y;
                ow.writeSignedInteger(ow.modal.placement_y);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0 when count_x == 1:
                ow.writeUnsignedInteger(3);
                ow.modal.y_dimension = count_y;
                ow.writeUnsignedInteger((uint)(count_y - 2));
                ow.modal.y_space = pitch.Y;
                ow.writeUnsignedInteger((uint)pitch.Y);
                break;
            case > 0 when count_y == 1:
                ow.writeUnsignedInteger(2);
                ow.modal.x_dimension = count_x;
                ow.writeUnsignedInteger((uint)(count_x - 2));
                ow.modal.x_space = pitch.X;
                ow.writeUnsignedInteger((uint)pitch.X);
                break;
            case > 0:
                ow.writeUnsignedInteger(1);
                ow.modal.x_dimension = count_x;
                ow.modal.y_dimension = count_y;
                ow.writeUnsignedInteger((uint)(count_x - 2));
                ow.writeUnsignedInteger((uint)(count_y - 2));
                ow.modal.x_space = pitch.X;
                ow.modal.y_space = pitch.Y;
                ow.writeUnsignedInteger((uint)pitch.X);
                ow.writeUnsignedInteger((uint)pitch.Y);
                break;
        }
    }

    public override bool isCellrefArray()
    {
        return pIsCellRefArray();
    }

    private bool pIsCellRefArray()
    {
        return true;
    }

    public override GCCell getCellref() {
        return pGetCellRef();
    }

    private GCCell pGetCellRef()
    {
        return cell_ref;
    }

    public override GeoLibPoint getPos()
    {
        return pGetPos();
    }

    private GeoLibPoint pGetPos()
    {
        return point;
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }


    private List<GCPolygon> pConvertToPolygons()
    {
        /*
        // This returns the contents of the reference.
        List<GCPolygon> tmp = new List<GCPolygon>();
        for (int element = 0; element < cell_ref.elementList.Count; element++)
        {
            List<GCPolygon> lp = cell_ref.elementList[element].convertToPolygons().ToList();
            // Above should work with nested arrays due to the recursion of the convertToPolygons override.
            tmp.AddRange(lp);
        }
        */

        List<GCPolygon> tmp = cell_ref.convertToPolygons();

        // Now we need to array the above.
        List<GCPolygon> ret = new();
        for (int x = 0; x < count_x; x++)
        {
            for (int y = 0; y < count_y; y++)
            {
                foreach (GCPolygon tp in tmp.Select(t => new GCPolygon(t)))
                {
                    tp.move(new GeoLibPoint(x * pitch.X, y * pitch.Y));
                    ret.Add(tp);
                }
            }
        }

#if !GCSINGLETHREADED
        Parallel.For(0, ret.Count, (poly, loopstate) =>
#else
            for (int poly = 0; poly < ret.Count; poly++)
#endif
            {
                ret[poly].rotate(trans.angle, point);
                ret[poly].scale(trans.mag);
            }
#if !GCSINGLETHREADED
        );
#endif

        return ret;
    }

}