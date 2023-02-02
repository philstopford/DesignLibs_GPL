using System;

namespace gds;

internal partial class gdsReader
{
    private void addCellRef()
    {
        cell_.addCellref();
        cell_.elementList[^1].setPos(new (modal.point_array[0].X, modal.point_array[0].Y));
        cell_.elementList[^1].setCellRef(drawing_.findCell(modal.sname));
        cell_.elementList[^1].setName(modal.sname);
        cell_.elementList[^1].rotate(modal.angle);
        cell_.elementList[^1].scale(modal.mag);
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].setMirrorx();
                break;
        }
    }

    private void addCellRefArray()
    {
        /*if (mirror_x){
            point=point_array.point(1);
            point.setY(-point.y());
            point_array.setPoint(1,point);
            point=point_array.point(2);
            point.setY(-point.y());
            point_array.setPoint(2,point);};*/
        cell_.addCellrefArray(drawing_.findCell(modal.sname), modal.point_array, modal.anzx, modal.anzy);
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].rotate(-modal.angle);
                break;
            default:
                cell_.elementList[^1].rotate(modal.angle);
                break;
        }
        cell_.elementList[^1].scale(modal.mag);
        cell_.elementList[^1].setName(modal.sname);
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].setMirrorx();
                break;
        }
    }

    private void addBox()
    {
        cell_.addBox(modal.point_array, modal.layer, modal.datatype);
    }

    private void addPolygon()
    {
        cell_.addPolygon(modal.point_array, modal.layer, modal.datatype);
    }

    private void addText()
    {
        modal.width = modal.width switch
        {
            1 => -10,
            _ => modal.width switch
            {
                0 => -10,
                _ => modal.mag switch
                {
                    <= 1 when modal.width == 0 => -10,
                    _ => modal.width
                }
            }
        };
        cell_.addText(modal.layer, modal.datatype, modal.point_array[0], modal.sname);
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].rotate(-modal.angle);
                break;
            default:
                cell_.elementList[^1].rotate(modal.angle);
                break;
        }
        cell_.elementList[^1].scale(modal.mag);
        cell_.elementList[^1].setName(modal.sname);
        switch (modal.mirror_x)
        {
            case true:
                cell_.elementList[^1].setMirrorx();
                break;
        }
        cell_.elementList[^1].setWidth(modal.width);
        cell_.elementList[^1].setPresentation(modal.presentation);
    }

    private void addPath()
    {
        switch (modal.point_array.Count)
        {
            case 1:
                cell_.addCircle(modal.layer, modal.datatype, modal.point_array[0], Convert.ToDouble(modal.width) / 2);
                break;
            default:
            {
                cell_.addPath(modal.point_array, modal.layer, modal.datatype);
                cell_.elementList[^1].setWidth(modal.width);
                if (modal.cap != 4)
                {
                    cell_.elementList[^1].setCap(modal.cap);
                }
                else
                {
                    cell_.elementList[^1].expandCaps(modal.beginExt, modal.endExt);
                }

                break;
            }
        }
    }
}