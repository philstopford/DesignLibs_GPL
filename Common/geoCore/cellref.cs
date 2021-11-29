using gds;
using geoLib;
using oasis;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace geoCoreLib;

public class GCCellref : GCElement
{
    public GeoLibPoint point { get; set; }
    private string name;
    public GCCell cell_ref { get; set; }
    public GCStrans trans { get; set; }

    public override void setName(string s)
    {
        pSetName(s);
    }

    private void pSetName(string s)
    {
        name = s;
    }

    public GCCellref(GCCell c, GeoLibPoint pos)
    {
        pGCCellref(c, pos);
    }

    private void pGCCellref(GCCell c, GeoLibPoint pos)
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

    private void pGCCellref()
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

    private void pSetPos(GeoLibPoint p)
    {
        point = new GeoLibPoint(p.X, p.Y);
    }

    public override GeoLibPoint getPos()
    {
        return pGetPos();
    }

    private GeoLibPoint pGetPos()
    {
        return point;
    }

    public override void minimum(GeoLibPoint p)
    {
        pMinimum(p);
    }

    private void pMinimum(GeoLibPoint p)
    {
        GeoLibPoint pos1 = new(p.X - point.X, p.Y - point.Y);
        pos1.Y = trans.mirror_x switch
        {
            true => -pos1.Y,
            _ => pos1.Y
        };
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
                break;
        }
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
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
                break;
        }
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos2.Y = -pos2.Y;
                pos2.X = -pos2.X;
                break;
        }
        pos1.Offset(point);
        pos2.Offset(point);
        p.X = Math.Min(p.X, pos2.X);
        p.X = Math.Min(p.X, pos1.X);
        p.Y = Math.Min(p.Y, pos2.Y);
        p.Y = Math.Min(p.Y, pos1.Y);
    }

    public override void scale(double factor)
    {
        pScale(factor);
    }

    private void pScale(double factor)
    {
        trans.mag = factor;
    }

    public override void rotate(double angle)
    {
        pRotate(angle);
    }

    private void pRotate(double angle)
    {
        trans.angle = angle;
    }

    public override void maximum(GeoLibPoint p)
    {
        pMaximum(p);
    }

    private void pMaximum(GeoLibPoint p)
    {
        GeoLibPoint pos1 = new(p.X - point.X, p.Y - point.Y);
        pos1.Y = trans.mirror_x switch
        {
            true => -pos1.Y,
            _ => pos1.Y
        };
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
                break;
        }
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
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos1.Y = -pos1.Y;
                pos1.X = -pos1.X;
                break;
        }
        switch (Math.Abs(trans.angle - 180))
        {
            case <= double.Epsilon:
                pos2.Y = -pos2.Y;
                pos2.X = -pos2.X;
                break;
        }
        pos1.Offset(point);
        pos2.Offset(point);
        p.X = Math.Max(p.X, pos2.X);
        p.X = Math.Max(p.X, pos1.X);
        p.Y = Math.Max(p.Y, pos2.Y);
        p.Y = Math.Max(p.Y, pos1.Y);
    }

    public override GCCell depend()
    {
        return pDepend();
    }

    private GCCell pDepend()
    {
        return cell_ref;
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
        point.Offset(p);
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
    }

    public override bool isCellref()
    {
        return pIsCellref();
    }

    private bool pIsCellref()
    {
        return true;
    }

    public override void setCellRef(GCCell cellRef)
    {
        pSetCellRef(cellRef);
    }

    private void pSetCellRef(GCCell cellRef)
    {
        cell_ref = cellRef;
    }

    public override GCCell getCellref()
    {
        return pGetCellRef();
    }

    private GCCell pGetCellRef()
    {
        return cell_ref;
    }

    public override void saveGDS(gdsWriter bw)
    {
        pSaveGDS(bw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        switch (cell_ref)
        {
            // Guard against incoming broken cellref information
            case null:
                return;
        }
        //SRef
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x0A);
        gw.bw.Write((byte)0);
        gw.writeString(cell_ref.cellName, 0x12);
        int strans_ = 0;
        switch (trans.mirror_x)
        {
            case true:
                strans_ |= 32768;
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
        //xy
        gw.bw.Write((ushort)(2 * 4 + 4));
        gw.bw.Write((byte)0x10);
        gw.bw.Write((byte)3);
        gw.bw.Write(point.X);
        gw.bw.Write(point.Y);
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write((byte)0x11);
        gw.bw.Write((byte)0);
    }

    public override void saveOASIS(oasWriter ow)
    {
        pSaveOASIS(ow);
    }

    private void pSaveOASIS(oasWriter ow)
    {
        ow.modal.absoluteMode = ow.modal.absoluteMode switch
        {
            false => true,
            _ => ow.modal.absoluteMode
        };
        byte info_byte = 128;  //explicid cellname;
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
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        /*
        List<GCPolygon> ret = new List<GCPolygon>();

        for (int element = 0; element < cell_ref.elementList.Count; element++)
        {
            List<GCPolygon> tmp = cell_ref.elementList[element].convertToPolygons();
            ret.AddRange(tmp);
        }
        */
        List<GCPolygon> ret = cell_ref.convertToPolygons();


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