using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoCoreLib;

public class GCDrawingfield
{
    public static double default_databaseunits = 1E-9;
    public static double default_userunits = 1E-3;

    public List<GCCell> cellList { get; set; }
    public int active_cell { get; set; }
    public double databaseunits { get; set; }
    public string libname { get; set; }
    public double userunits { get; set; }

    public short modyear { get; set; }
    public short modmonth { get; set; }
    public short modday { get; set; }
    public short modhour { get; set; }
    public short modmin { get; set; }
    public short modsec { get; set; }

    public short accyear { get; set; }
    public short accmonth { get; set; }
    public short accday { get; set; }
    public short acchour { get; set; }
    public short accmin { get; set; }
    public short accsec { get; set; }

    public GCDrawingfield copy()
    {
        return (GCDrawingfield)MemberwiseClone();
    }

    public void reset()
    {
        pReset();
    }

    private void pReset()
    {
        cellList = new List<GCCell>();
        active_cell = 0;
        databaseunits = default_databaseunits;
        userunits = default_userunits;
        libname = "noname";
        DateTime now = DateTime.Now;
        modyear = (short)now.Year;
        modmonth = (short)now.Month;
        modday = (short)now.Day;
        modhour = (short)now.Hour;
        modmin = (short)now.Minute;
        modsec = (short)now.Second;
        accyear = (short)now.Year;
        accmonth = (short)now.Month;
        accday = (short)now.Day;
        acchour = (short)now.Hour;
        accmin = (short)now.Minute;
        accsec = (short)now.Second;
    }

    public GCDrawingfield(string name)
    {
        pGCDrawingfield(name);
    }

    private void pGCDrawingfield(string name)
    {
        pReset();
        libname = name;
    }

    public GCDrawingfield(GCDrawingfield drawing)
    {
        pGCDrawingfield(drawing);
    }

    private void pGCDrawingfield(GCDrawingfield drawing)
    {
        if (drawing == null)
        {
            reset();
            return;
        }
        try
        {
            cellList = drawing.cellList.ToList();
        }
        catch
        {
            cellList = new();
        }

        accyear = drawing.accyear;
        accmonth = drawing.accmonth;
        accday = drawing.accday;
        acchour = drawing.acchour;
        accmin = drawing.accmin;
        accsec = drawing.accsec;
        modyear = drawing.modyear;
        modmonth = drawing.modmonth;
        modday = drawing.modday;
        modhour = drawing.modhour;
        modmin = drawing.modmin;
        modsec = drawing.modsec;
        active_cell = drawing.active_cell;
        databaseunits = drawing.databaseunits;
        libname = drawing.libname;
        userunits = drawing.userunits;
    }
    
    public GCCell findCell(string s)
    {
        return pFindCell(s);
    }

    public int findCellIndex(string s)
    {
        return pFindCell_serial(s);
    }

    private int pFindCell_serial(string s)
    {
        int found = -1;
        for (int i = 0; i < cellList.Count; i++)
        {
            string name = cellList[i].cellName;
            if (name != s)
            {
                continue;
            }

            found = i;
            break;
        }

        return found;
    }

    private GCCell pFindCell(string s)
    {
        int found = pFindCell_serial(s);
        return found > -1 ? cellList[found] : null;
    }

    public GCCell findCell_par(string s)
    {
        return pFindCell_par(s);
    }

    private GCCell pFindCell_par(string s)
    {
        int found = -1;
        for (int i = 0; i < cellList.Count; i++)
        {
            string name = cellList[i].cellName;
            if (name != s)
            {
                continue;
            }

            found = i;
            break;
        }
        return found > -1 ? cellList[found] : null;
    }

    public void addCells(int count)
    {
        pAddCells(count);
    }

    private void pAddCells(int count)
    {
        cellList = cellList switch
        {
            null => new List<GCCell>(),
            _ => cellList
        };

        cellList.AddRange(new GCCell[count].ToList());
    }

    public GCCell addCell()
    {
        return pAddCell();
    }

    private GCCell pAddCell()
    {
        cellList = cellList switch
        {
            null => new List<GCCell>(),
            _ => cellList
        };
        cellList.Add(new GCCell());
        return cellList[^1];
    }

    public GCCell findTopCell()
    {
        return pFindTopCell();
    }

    private GCCell pFindTopCell()
    {
        List<GCCell> cell_list = cellList.ToList();
        List<GCCell> celllist = cellList.ToList();

        int i = 0;
        int j = 0;
        for (i = 0; i < cell_list.Count; i++)
        {
            bool topcell = true;
            for (j = 0; j < celllist.Count; j++)
            {
                if (cell_list[i] == celllist[i])
                {
                    continue;
                }

                if (celllist[i].depend(cell_list[i]))
                {
                    topcell = false;
                }
            }
            switch (topcell)
            {
                case true:
                    return cell_list[i];
            }
        }
        return cellList[j];
    }

    public GCCell findTopCell_par()
    {
        return pFindTopCell_par();
    }

    private GCCell pFindTopCell_par()
    {
        List<GCCell> cell_list = cellList.ToList();
        List<GCCell> celllist = cellList.ToList();
        int i;
        int j = 0;
        for (i = 0; i < cell_list.Count; i++)
        {
            bool topcell = true;
            for (j = 0; j < celllist.Count; j++)
            {
                if (cell_list[i] == celllist[i])
                {
                    continue;
                }

                if (celllist[i].depend(cell_list[i]))
                {
                    topcell = false;
                }
            }
            switch (topcell)
            {
                case true:
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

    private void pResize(double scale)
    {
        int cellListCount = cellList.Count;
#if !GCSINGLETHREADED
        Parallel.For(0, cellListCount, i =>
#else
            for (int i = 0; i < cellListCount; i++)
#endif
            {
                if (cellList[i].elementList != null)
                {
                    cellList[i].resize(scale);
                }
            }
#if !GCSINGLETHREADED
        );
#endif

        databaseunits /= scale;
        userunits *= scale;
    }

    public List<List<GCPolygon>> convertToPolygons(int layer = -1, int datatype = -1, List<string> cells = null)
    {
        return pConvertToPolygons(layer, datatype, cells);
    }

    private List<List<GCPolygon>> pConvertToPolygons(int layer = -1, int datatype = -1, List<string> cells = null)
    {
        List<GCPolygon>[] ret = Array.Empty<List<GCPolygon>>();

        if (cells == null)
        {
            int cellCount = cellList.Count;
            ret = new List<GCPolygon>[cellCount];

#if !GCSINGLETHREADED
            ParallelOptions po = new();
            Parallel.For(0, cellCount, po, i =>
#else
            for (int i = 0; i < cellList.Count; i++)
#endif
                {
                    ret[i] = cellList[i].convertToPolygons(layer: layer, datatype: datatype);
                }
#if !GCSINGLETHREADED
            );
#endif
        }
        else
        {
            ret = new List<GCPolygon>[cells.Count];
#if !GCSINGLETHREADED
            ParallelOptions po = new();
            Parallel.For(0, ret.Length, po, c =>
#else
            for (int c = 0; c < ret.Length; c++)
#endif
            {
                string cell = cells[c];
                int cell_index = pFindCell_serial(cell);
                if (cell_index != -1)
                {
                    ret[c] = cellList[cell_index].convertToPolygons(layer: layer, datatype: datatype);
                }
                else
                {
                    ret[c] = new List<GCPolygon>();
                }
            }
#if !GCSINGLETHREADED
            );
#endif
        }

        // Apply scaling...if needed.
        if (databaseunits != default_databaseunits)
        {
            double grid_scaling = databaseunits / default_databaseunits;
#if !GCSINGLETHREADED
            ParallelOptions po1 = new();
            Parallel.For(0, ret.Length, po1, cindex =>
#else
            for (int cindex = 0; cindex < ret.Length; cindex++)
#endif
            {
#if !GCSINGLETHREADED
                Parallel.For(0, ret[cindex].Count, po1, i =>
#else
                for (int i = 0; i < ret[cindex].Count; i++)
#endif
                {
                    ret[cindex][i].resize(grid_scaling);
                }
#if !GCSINGLETHREADED
                );
#endif
            }
#if !GCSINGLETHREADED
            );
#endif
        }

        return ret.ToList();
    }

}