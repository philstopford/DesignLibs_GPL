using System;

namespace geoCoreLib
{
    [Serializable]
    public class GeoCoreHandler
    {
        Boolean changed;

        public bool isChanged()
        {
            return pIsChanged();
        }

        bool pIsChanged()
        {
            return changed;
        }

        public void setChanged(bool value)
        {
            pSetChanged(value);
        }

        void pSetChanged(bool value)
        {
            changed = value;
        }

        Boolean valid;
        public string errormsg;
        public bool isValid()
        {
            return pIsValid();
        }

        bool pIsValid()
        {
            return valid;
        }

        public void setValid(bool value)
        {
            pSetValid(value);
        }

        void pSetValid(bool value)
        {
            valid = value;
        }

        string filename;
        public string getFilename()
        {
            return pGetFilename();
        }

        string pGetFilename()
        {
            return filename;
        }

        public void setFilename(string name)
        {
            pSetFilename(name);
        }

        void pSetFilename(string name)
        {
            filename = name;
        }

        GeoCore geo;
        public GeoCore getGeo()
        {
            return pGetCore();
        }

        GeoCore pGetCore()
        {
            return geo;
        }

        public GeoCoreHandler()
        {
            init();
        }

        void init()
        {
            // stubs
            changed = false;
            filename = "";
            geo = new GeoCore();
            valid = geo.isValid();
            errormsg = geo.errormsg;
        }

        public void reset()
        {
            pReset();
        }

        void pReset()
        {
            changed = true;
            geo.reset();
            valid = geo.isValid();
            errormsg = geo.errormsg;
            filename = "";
        }

        public void readValues(GeoCoreHandler sourceGeoCoreHandler)
        {
            pReadValues(sourceGeoCoreHandler);
        }

        void pReadValues(GeoCoreHandler sourceGeoCoreHandler)
        {
            valid = sourceGeoCoreHandler.valid;
            errormsg = sourceGeoCoreHandler.errormsg;
            filename = sourceGeoCoreHandler.filename;
            if (geo == null)
            {
                geo = new GeoCore();
            }
            geo.readValues(sourceGeoCoreHandler.geo);
        }

        public void updateGeoCoreHandler(string filename, GeoCore.fileType type)
        {
            pUpdateGeoCoreHandler(filename, type);
        }

        void pUpdateGeoCoreHandler(string filename, GeoCore.fileType type)
        {
            valid = true;
            this.filename = filename;

            geo.updateGeoCore(filename, type);
            valid = geo.isValid();
            errormsg = geo.errormsg;
        }

        public void setPoints(Int32 structureIndex, Int32 ldIndex)
        {
            pSetPoints(structureIndex, ldIndex);
        }

        void pSetPoints(Int32 structureIndex, Int32 ldIndex)
        {
            if (changed)
            {
                geo.updateGeometry(structureIndex, ldIndex);
                changed = false;
            }
        }
    }
}
