using gds;
using geoWrangler;
using oasis;
using System;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoCoreLib;

public class GCPath : GCElement
{
    public Path64 pointarray;
    private int width;
    private int cap;
    private bool isRound;

    public GCPath()
    {
        pGCPath();
    }

    private void pGCPath()
    {
        pointarray = new();
    }

    public GCPath(Path64 points, int layer, int datatype)
    {
        pGCPath(points, layer, datatype);
    }

    private void pGCPath(Path64 points, int layer, int datatype)
    {
        layer_nr = layer;
        datatype_nr = datatype;
        pointarray = new(points);
        width = 0;
        cap = 0;
        clean();
    }

    public override bool correct()
    {
        return pCorrect();
    }

    private bool pCorrect()
    {
        switch (pointarray.Count)
        {
            case < 2:
            case 2 when pointarray[0] == pointarray[1] && cap == 0:
                return false;
            default:
                return true;
        }
    }

    public override void minimum(ref Point64 pos)
    {
        pMinimum(ref pos);
    }

    private void pMinimum(ref Point64 pos)
    {
        Point64 p = GeoWrangler.getMinimumPoint(pointarray);
        pos.X = Math.Min(pos.X, p.X);
        pos.Y = Math.Min(pos.Y, p.Y);
    }

    public override void maximum(ref Point64 pos)
    {
        pMaximum(ref pos);
    }

    private void pMaximum(ref Point64 pos)
    {
        Point64 p = GeoWrangler.getMaximumPoint(pointarray);
        pos.X = Math.Max(pos.X, p.X);
        pos.Y = Math.Max(pos.Y, p.Y);
    }

    public override void moveSelect(Point64 pos)
    {
        pMoveSelect(pos);
    }

    private void pMoveSelect(Point64 pos)
    {
        switch (select)
        {
            case true:
            {
                pointarray = GeoWrangler.move(pointarray, pos.X, pos.Y);
                break;
            }
        }
    }

    public override void move(Point64 pos)
    {
        pMove(pos);
    }

    private void pMove(Point64 pos)
    {
        pointarray = GeoWrangler.move(pointarray, pos.X, pos.Y);
    }

    public override void resize(double size)
    {
        pResize(size);
    }

    private void pResize(double size)
    {
        width = (int)(size * width);
        pointarray = GeoWrangler.resize(pointarray, size);
    }

    public override void clean()
    {
        pClean();
    }

    private void pClean()
    {
        for (int i = 0; i < pointarray.Count - 1; i++)
        {
            switch (pointarray.Count)
            {
                case <= 2:
                    return; //no area
            }

            if (pointarray[i] != pointarray[i + 1])
            {
                continue;
            }

            deletePoint(i + 1);
            i--;
        }
    }

    public void deletePoint(int pos)
    {
        pDeletePoint(pos);
    }

    private void pDeletePoint(int pos)
    {
        pointarray.RemoveAt(pos);
    }

    public override void setWidth(int width_)
    {
        pSetWidth(width_);
    }

    private void pSetWidth(int width_)
    {
        width = width_;
    }

    public override void setCap(int cap_)
    {
        pSetCap(cap_);
    }

    private void pSetCap(int cap_)
    {
        cap = cap_;
    }

    public override Path64 getPath()
    {
        return pGetPath();
    }

    private Path64 pGetPath()
    {
        return pointarray;
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        List<GCPolygon> ret = new();
        Path64 tmp = GeoWrangler.inflatePath(pointarray, (double)width / 2);
        ret.Add(new GCPolygon(tmp, layer_nr, datatype_nr));
        return ret;
    }

    public void setRound(bool val)
    {
        pSetRound(val);
    }

    private void pSetRound(bool val)
    {
        isRound = val;
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        // path 
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sPATH);
        //layer
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sLAYER);
        gw.bw.Write((short)layer_nr);
        //datatype
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sDATATYPE);
        gw.bw.Write((short)datatype_nr);
        //pathtype
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sPATHTYPE);
        switch (isRound)
        {
            case false:
                gw.bw.Write((short)cap);
                break;
            default:
                gw.bw.Write((short)1);
                break;
        }
        //width
        gw.bw.Write((ushort)8);
        gw.bw.Write(gdsValues.sWIDTH);
        gw.bw.Write(width);

        int i = pointarray.Count;
        i = i switch
        {
            > 8191 => 8191,
            _ => i
        };
        //xy 
        int val = i * 2 * 4 + 4;
        gw.bw.Write((ushort)val);
        gw.bw.Write(gdsValues.sXY);
        for (int k = 0; k < i; k++)
        {
            gw.bw.Write((int)pointarray[k].X);
            gw.bw.Write((int)pointarray[k].Y);
        }
        // endel
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sENDEL);
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
        if (width / 2 != ow.modal.geometry_w)
        {
            info_byte += 64;
        }

        ow.writeUnsignedInteger(22);
        ow.writeRaw(info_byte);
        switch (info_byte & 1)
        {
            case > 0:
                ow.modal.layer = layer_nr;
                ow.writeUnsignedInteger((uint)ow.modal.layer);
                break;
        }
        switch (info_byte & 2)
        {
            case > 0:
                ow.modal.datatype = datatype_nr;
                ow.writeUnsignedInteger((uint)datatype_nr);
                break;
        }
        switch (info_byte & 64)
        {
            case > 0:
            {
                ow.modal.geometry_w = ow.modal.geometry_w switch
                {
                    < 0 => -ow.modal.geometry_w,
                    // Yes, indeed. Oasis is limited to path widths that are factors of 2. Strange, but there you are.
                    _ => width / 2
                };
                ow.writeUnsignedInteger((uint)ow.modal.geometry_w);
                break;
            }
        }

        switch (cap)
        {
            // 2 means the end gets extended by the half the path width.
            case 0:
                ow.writeUnsignedInteger(5);
                break;
            case 2:
                ow.writeUnsignedInteger(10);
                break;
            default:
                ow.writeUnsignedInteger(5);
                // Round caps are not possible in OASIS files.;
                break;
        }
        ow.writePointArray(pointarray, false);
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.geometry_x = (int)pointarray[0].X;
                ow.writeSignedInteger(ow.modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                ow.modal.geometry_y = (int)pointarray[0].Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
                break;
        }
    }

    public override bool isPath()
    {
        return true;
    }
}