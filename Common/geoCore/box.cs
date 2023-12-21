using System;
using gds;
using geoLib;
using oasis;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoCoreLib;

public class GCBox : GCElement
{
    private GeoLibRectangle rect;

    public GCBox()
    {
    }

    public GCBox(GeoLibRectangle r, int layer, int datatype)
    {
        pGCBox(r, layer, datatype);
    }

    private void pGCBox(GeoLibRectangle r, int layer, int datatype)
    {
        rect = r;
        layer_nr = layer;
        datatype_nr = datatype;
    }

    public GCBox(int x_, int y_, int b_, int h_, int layer, int datatype)
    {
        pGCBox(x_, y_, b_, h_, layer, datatype);
    }

    private void pGCBox(int x_, int y_, int b_, int h_, int layer, int datatype)
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

    private bool pCorrect()
    {
        if (rect.Width == 0)
        {
            return false;
        }
        return rect.Height != 0;
    }

    public override void maximum(ref Point64 p)
    {
        pMaximum(ref p);
    }

    private void pMaximum(ref Point64 p)
    {
        p.X = rect.Location.X + rect.Width;
        p.Y = rect.Location.Y + rect.Height;
    }

    public override void minimum(ref Point64 p)
    {
        pMinimum(ref p);
    }

    private void pMinimum(ref Point64 p)
    {
        p.X = rect.Location.X;
        p.Y = rect.Location.Y;
    }

    public override void move(Point64 pos)
    {
        pMove(pos);
    }

    private void pMove(Point64 pos)
    {
        rect.Offset(pos);
    }

    public override void moveSelect(Point64 p)
    {
        pMoveSelect(p);
    }

    private void pMoveSelect(Point64 p)
    {
        switch (select)
        {
            case true:
                rect.Offset(p);
                break;
        }
    }

    public override void resize(double size)
    {
        pResize(size);
    }

    private void pResize(double size)
    {
        rect = new GeoLibRectangle(
            (int)(rect.Location.X * size),
            (int)(rect.Location.Y * size),
            (int)(rect.Width * size),
            (int)(rect.Height * size)
        );
    }

    public override List<GCPolygon> convertToPolygons()
    {
        return pConvertToPolygons();
    }

    private List<GCPolygon> pConvertToPolygons()
    {
        List<GCPolygon> ret = new();
        Path64 points = Helper.initedPath64(5);
        Int64 left = rect.Location.X;
        Int64 bottom = rect.Location.Y;
        Int64 right = left + rect.Width;
        Int64 top = bottom + rect.Height;
        points[0] = new (left, top);
        points[1] = new (right, top);
        points[2] = new (right, bottom);
        points[3] = new (left, bottom);
        points[4] = new (left, top);
        ret.Add(new GCPolygon(points, layer_nr, datatype_nr));
        return ret;
    }

    public override void saveGDS(gdsWriter gw)
    {
        pSaveGDS(gw);
    }

    private void pSaveGDS(gdsWriter gw)
    {
        int left = (int)rect.Location.X;
        int bottom = (int)rect.Location.Y;
        int right = left + rect.Width;
        int top = bottom + rect.Height;

        // box
        gw.bw.Write((ushort)4);
        gw.bw.Write(gdsValues.sBOUNDARY);
        //layer
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sLAYER);
        gw.bw.Write((short)layer_nr);
        //boxtype
        gw.bw.Write((ushort)6);
        gw.bw.Write(gdsValues.sDATATYPE);
        gw.bw.Write((short)datatype_nr);
        //xy 
        // 5 points (last must match first). 2 values per point. 4 bytes per value.
        gw.bw.Write((ushort)(4 + 5 * 2 * 4));
        gw.bw.Write(gdsValues.sXY);
        gw.bw.Write(left);
        gw.bw.Write(top);
        gw.bw.Write(right);
        gw.bw.Write(top);
        gw.bw.Write(right);
        gw.bw.Write(bottom);
        gw.bw.Write(left);
        gw.bw.Write(bottom);
        gw.bw.Write(left);
        gw.bw.Write(top);
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

        byte info_byte = 0;
        if (layer_nr != ow.modal.layer)
        {
            info_byte += 1;
        }
        if (datatype_nr != ow.modal.datatype)
        {
            info_byte += 2;
        }
        if (rect.Location.X != ow.modal.geometry_x)
        {
            info_byte += 16;
        }
        if (rect.Location.Y != ow.modal.geometry_y)
        {
            info_byte += 8;
        }
        if (rect.Height != ow.modal.geometry_h)
        {
            info_byte += 32;
        }
        if (rect.Width != ow.modal.geometry_w)
        {
            info_byte += 64;
        }
        if (rect.Width == rect.Height)
        {
            info_byte += 128;
            info_byte = (byte)(info_byte & (255 - 32));
        }
        ow.writeUnsignedInteger(20);
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
                ow.modal.geometry_w = rect.Width;
                ow.writeUnsignedInteger((uint)ow.modal.geometry_w);

                break;
            }
        }

        switch (info_byte & 128)
        {
            case > 0:
                ow.modal.geometry_h = ow.modal.geometry_w;
                break;
            default:
                ow.modal.geometry_h = rect.Height;
                break;
        }

        switch (info_byte & 32)
        {
            case > 0:
                ow.writeUnsignedInteger((uint)ow.modal.geometry_h);
                break;
        }
        switch (info_byte & 16)
        {
            case > 0:
                ow.modal.geometry_x = (int)rect.Location.X;
                ow.writeSignedInteger(ow.modal.geometry_x);
                break;
        }
        switch (info_byte & 8)
        {
            case > 0:
                ow.modal.geometry_y = (int)rect.Location.Y;
                ow.writeSignedInteger(ow.modal.geometry_y);
                break;
        }
    }

    public override bool isBox()
    {
        return true;
    }

    public override Point64 getPos()
    {
        return new (rect.Location);
    }

    public override int getWidth()
    {
        return rect.Width;
    }
    
    public override int getHeight()
    {
        return rect.Height;
    }
}