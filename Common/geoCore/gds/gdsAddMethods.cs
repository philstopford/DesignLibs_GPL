using System;
using Clipper2Lib;

namespace gds;

internal partial class gdsReader
{
    private void addCellRef()
    {
        cell_.addCellref();
        cell_.elementList[^1].setPos(new Point64(modal.point_array[0].X, modal.point_array[0].Y));
        cell_.elementList[^1].setCellRef(drawing_.findCell(modal.sname));
        cell_.elementList[^1].setName(modal.sname);
        cell_.elementList[^1].rotate(modal.angle);
        cell_.elementList[^1].scale(modal.mag);
        if (modal.mirror_x)
        {
            cell_.elementList[^1].setMirrorx();
        }
    }

    private void addCellRefArray()
    {
        (modal.point_array[1], modal.point_array[2]) = (modal.point_array[2], modal.point_array[1]);
        cell_.addCellrefArray(drawing_.findCell(modal.sname), modal.point_array, modal.anzx, modal.anzy);
        cell_.elementList[^1].rotate(modal.angle);
        cell_.elementList[^1].scale(modal.mag);
        cell_.elementList[^1].setName(modal.sname);
        if (modal.mirror_x)
        {
            cell_.elementList[^1].setMirrorx();
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
        switch (modal.width)
        {
            case 1:
            case 0:
                modal.width = -10;
                break;
            default:
            {
                if (modal is { mag: <= 1, width: 0 })
                {
                    modal.width = -10;
                }
                else
                {
                    modal.width = modal.width;
                }

                break;
            }
        }

        cell_.addText(modal.layer, modal.datatype, modal.point_array[0], modal.sname);
        cell_.elementList[^1].rotate(modal.angle);
        cell_.elementList[^1].scale(modal.mag);
        cell_.elementList[^1].setName(modal.sname);
        if (modal.mirror_x)
        {
            cell_.elementList[^1].setMirrorx();
        }

        cell_.elementList[^1].setWidth(modal.width);
        cell_.elementList[^1].setPresentation(modal.presentation);
    }

    private void addPath()
    {
        if (modal.point_array.Count == 1)
        {
            cell_.addCircle(modal.layer, modal.datatype, modal.point_array[0], Convert.ToDouble(modal.width) / 2);
        }
        else
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
        }
    }
}