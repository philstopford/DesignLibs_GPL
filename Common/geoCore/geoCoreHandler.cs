using System;
using System.Collections.Generic;

namespace geoCoreLib;

[Serializable]
public class GeoCoreHandler
{
    private bool changed;

    public bool isChanged()
    {
        return pIsChanged();
    }

    private bool pIsChanged()
    {
        return changed;
    }

    public void setChanged(bool value)
    {
        pSetChanged(value);
    }

    private void pSetChanged(bool value)
    {
        changed = value;
    }

    private bool valid;
    public List<string> error_msgs;
    public bool isValid()
    {
        return pIsValid();
    }

    private bool pIsValid()
    {
        return valid;
    }

    public void setValid(bool value)
    {
        pSetValid(value);
    }

    private void pSetValid(bool value)
    {
        valid = value;
    }

    private string filename;
    public string getFilename()
    {
        return pGetFilename();
    }

    private string pGetFilename()
    {
        return filename;
    }

    public void setFilename(string name)
    {
        pSetFilename(name);
    }

    private void pSetFilename(string name)
    {
        filename = name;
    }

    private GeoCore geo;
    public GeoCore getGeo()
    {
        return pGetCore();
    }

    private GeoCore pGetCore()
    {
        return geo;
    }

    public GeoCoreHandler()
    {
        init();
    }

    private void init()
    {
        // stubs
        changed = false;
        filename = "";
        geo = new GeoCore();
        valid = geo.isValid();
        error_msgs = geo.error_msgs;
    }

    public void reset()
    {
        pReset();
    }

    private void pReset()
    {
        changed = true;
        geo.reset();
        valid = geo.isValid();
        error_msgs = geo.error_msgs;
        filename = "";
    }

    public void readValues(GeoCoreHandler sourceGeoCoreHandler)
    {
        pReadValues(sourceGeoCoreHandler);
    }

    private void pReadValues(GeoCoreHandler sourceGeoCoreHandler)
    {
        valid = sourceGeoCoreHandler.valid;
        error_msgs = sourceGeoCoreHandler.error_msgs;
        filename = sourceGeoCoreHandler.filename;
        geo = geo switch
        {
            null => new GeoCore(),
            _ => geo
        };
        geo.readValues(sourceGeoCoreHandler.geo);
    }

    public void updateGeoCoreHandler(string filename_, GeoCore.fileType type)
    {
        pUpdateGeoCoreHandler(filename_, type);
    }

    private void pUpdateGeoCoreHandler(string filename_, GeoCore.fileType type)
    {
        valid = true;
        filename = filename_;

        geo.updateGeoCore(filename, type);
        valid = geo.isValid();
        error_msgs = geo.error_msgs;
    }

    public void setPoints(int structureIndex, int ldIndex)
    {
        pSetPoints(structureIndex, ldIndex);
    }

    private void pSetPoints(int structureIndex, int ldIndex)
    {
        switch (changed)
        {
            case true:
                geo.updateGeometry(structureIndex, ldIndex);
                changed = false;
                break;
        }
    }
}