using System;
using System.Collections.Generic;
using System.Linq;

namespace geoCoreLib
{
    public class GCDrawingfield
    {
        public List<GCCell> cellList { get; set; }
        public int active_cell { get; set; }
        public double databaseunits { get; set; }
        public string libname { get; set; }
        public double userunits { get; set; }

        public Int16 modyear { get; set; }
        public Int16 modmonth { get; set; }
        public Int16 modday { get; set; }
        public Int16 modhour { get; set; }
        public Int16 modmin { get; set; }
        public Int16 modsec { get; set; }

        public Int16 accyear { get; set; }
        public Int16 accmonth { get; set; }
        public Int16 accday { get; set; }
        public Int16 acchour { get; set; }
        public Int16 accmin { get; set; }
        public Int16 accsec { get; set; }

        public GCDrawingfield copy()
        {
            return (GCDrawingfield)MemberwiseClone();
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            cellList = new List<GCCell>();
            active_cell = 0;
            databaseunits = 1E-9;
            userunits = 1E-3;
            libname = "noname";
            modyear = 0;
            modmonth = 0;
            modday = 0;
            modhour = 0;
            modmin = 0;
            modsec = 0;
            accyear = 0;
            accmonth = 0;
            accday = 0;
            acchour = 0;
            accmin = 0;
            accsec = 0;
        }

        public GCDrawingfield(string name)
        {
            pGCDrawingfield(name);
        }

        void pGCDrawingfield(string name)
        {
            pReset();
            libname = name;
        }

        public GCCell findCell(string s)
        {
            return pFindCell(s);
        }

        public int findCellIndex(string s)
        {
            return pFindCell_serial(s);
        }

        int pFindCell_serial(string s)
        {
            int found = -1;
            for (int i = 0; i < cellList.Count; i++)
            {
                string name = cellList[i].cellName;
                if (name == s)
                {
                    found = i;
                    break;
                }
            }

            return found;
        }

        GCCell pFindCell(string s)
        {
            int found = pFindCell_serial(s);
            return found > -1 ? cellList[found] : null;
        }

        public GCCell findCell_par(string s)
        {
            return pFindCell_par(s);
        }

        GCCell pFindCell_par(string s)
        {
            int found = -1;
            for (int i = 0; i < cellList.Count; i++)
            {
                string name = cellList[i].cellName;
                if (name == s)
                {
                    found = i;
                    break;
                }
            }
            return found > -1 ? cellList[found] : null;
        }

        public void addCells(int count)
        {
            pAddCells(count);
        }

        void pAddCells(int count)
        {
            if (cellList == null)
            {
                cellList = new List<GCCell>();
            }
            cellList.AddRange(new GCCell[count].ToList());
        }

        public GCCell addCell()
        {
            return pAddCell();
        }

        GCCell pAddCell()
        {
            if (cellList == null)
            {
                cellList = new List<GCCell>();
            }
            cellList.Add(new GCCell());
            return cellList[^1];
        }

        public GCCell findTopCell()
        {
            return pFindTopCell();
        }

        GCCell pFindTopCell()
        {
            List<GCCell> cell_list = cellList.ToList();
            List<GCCell> celllist = cellList.ToList();

            int i = 0;
            int j = 0;
            for (i = 0; i < cell_list.Count(); i++)
            {
                bool topcell = true;
                for (j = 0; j < celllist.Count; j++)
                {
                    if (cell_list[i] != celllist[i])
                    {
                        if (celllist[i].depend(cell_list[i]))
                        {
                            topcell = false;
                        }
                    }
                }
                if (topcell)
                {
                    return cell_list[i];
                }
            }
            return cellList[j];
        }

        public GCCell findTopCell_par()
        {
            return pFindTopCell_par();
        }

        GCCell pFindTopCell_par()
        {
            List<GCCell> cell_list = cellList.ToList();
            List<GCCell> celllist = cellList.ToList();
            int i = 0;
            int j = 0;
            for (i = 0; i < cell_list.Count(); i++)
            {
                bool topcell = true;
                for (j = 0; j < celllist.Count; j++)
                {
                    if (cell_list[i] != celllist[i])
                    {
                        if (celllist[i].depend(cell_list[i]))
                        {
                            topcell = false;
                        }
                    }
                }
                if (topcell)
                {
                    return cell_list[i];
                }
            }
            return cellList[j];
        }

        // Experimental, to permit adapting to unit changes.
        public void resize(double scale)
        {
            pResize(scale);
        }

        void pResize(double scale)
        {
            int cellListCount = cellList.Count;
#if GCTHREADED
            Parallel.For(0, cellListCount, (i) =>
#else
            for (int i = 0; i < cellListCount; i++)
#endif
            {
                if (cellList[i].elementList != null)
                {
                    cellList[i].resize(scale);
                }
            }
#if GCTHREADED
            );
#endif

            databaseunits /= scale;
            userunits *= scale;
        }

        public List<List<GCPolygon>> convertToPolygons(int layer = -1, int datatype = -1)
        {
            return pConvertToPolygons(layer, datatype);
        }

        List<List<GCPolygon>> pConvertToPolygons(int layer = -1, int datatype = -1)
        {
            int cellCount = cellList.Count;
            List<GCPolygon>[] ret = new List<GCPolygon>[cellCount];

#if GCTHREADED
            ParallelOptions po = new ParallelOptions();
            Parallel.For(0, cellCount, po, i =>
#else
            for (int i = 0; i < cellList.Count; i++)
#endif
            {
                List<GCPolygon> cellPolys = cellList[i].convertToPolygons(layer:layer, datatype:datatype);
                ret[i] = cellPolys.ToList();
            }
#if GCTHREADED
            );
#endif
            return ret.ToList();
        }

    }
}
