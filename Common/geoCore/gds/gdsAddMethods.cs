using geoLib;

namespace gds
{
    partial class gdsReader
    {
        void addCellRef()
        {
            cell_.addCellref();
            cell_.elementList[cell_.elementList.Count - 1].setPos(new GeoLibPoint(modal.point_array[0].X, modal.point_array[0].Y));
            cell_.elementList[cell_.elementList.Count - 1].setCellRef(drawing_.findCell(modal.sname));
            cell_.elementList[cell_.elementList.Count - 1].setName(modal.sname);
            cell_.elementList[cell_.elementList.Count - 1].rotate(modal.angle);
            cell_.elementList[cell_.elementList.Count - 1].scale(modal.mag);
            if (modal.mirror_x)
            {
                cell_.elementList[cell_.elementList.Count - 1].setMirrorx();
            }
        }

        void addCellRefArray()
        {
            /*if (mirror_x){
				point=point_array.point(1);
				point.setY(-point.y());
				point_array.setPoint(1,point);
				point=point_array.point(2);
				point.setY(-point.y());
				point_array.setPoint(2,point);};*/
            cell_.addCellrefArray(drawing_.findCell(modal.sname), modal.point_array, modal.anzx, modal.anzy);
            if (modal.mirror_x)
            {
                cell_.elementList[cell_.elementList.Count - 1].rotate(-modal.angle);
            }
            else
            {
                cell_.elementList[cell_.elementList.Count - 1].rotate(modal.angle);
            }
            cell_.elementList[cell_.elementList.Count - 1].scale(modal.mag);
            cell_.elementList[cell_.elementList.Count - 1].setName(modal.sname);
            if (modal.mirror_x)
            {
                cell_.elementList[cell_.elementList.Count - 1].setMirrorx();
            }
        }

        void addBox()
        {
            cell_.addBox(modal.point_array, modal.layer, modal.datatype);
        }

        void addPolygon()
        {
            cell_.addPolygon(modal.point_array, modal.layer, modal.datatype);
        }

        void addText()
        {
            if (modal.mag <= 1 && modal.width == 0)
            {
                modal.width = -10;
            }
            if (modal.width == 0)
            {
                modal.width = -10;
            }
            if (modal.width == 1)
            {
                modal.width = -10;
            }
            cell_.addText(modal.layer, modal.datatype, modal.point_array[0], modal.sname);
            if (modal.mirror_x)
            {
                cell_.elementList[cell_.elementList.Count - 1].rotate(-modal.angle);
            }
            else
            {
                cell_.elementList[cell_.elementList.Count - 1].rotate(modal.angle);
            }
            cell_.elementList[cell_.elementList.Count - 1].scale(modal.mag);
            cell_.elementList[cell_.elementList.Count - 1].setName(modal.sname);
            if (modal.mirror_x)
            {
                cell_.elementList[cell_.elementList.Count - 1].setMirrorx();
            }
            cell_.elementList[cell_.elementList.Count - 1].setWidth(modal.width);
            cell_.elementList[cell_.elementList.Count - 1].setPresentation(modal.presentation);
        }

        void addPath()
        {
            if (modal.point_array.Length == 1)
            {
                cell_.addCircle(modal.layer, modal.datatype, modal.point_array[0], modal.width / 2);
            }
            else
            {
                cell_.addPath(modal.point_array, modal.layer, modal.datatype);
                cell_.elementList[cell_.elementList.Count - 1].setWidth(modal.width);
                if (modal.cap != 4)
                {
                    cell_.elementList[cell_.elementList.Count - 1].setCap(modal.cap);
                }
                else
                {
                    cell_.elementList[cell_.elementList.Count - 1].expandCaps(modal.beginExt, modal.endExt);
                }
            }
        }
    }
}
