﻿using System;
using System.Collections.Generic;
using entropyRNG;
using utility;

namespace LHC
{
    public class LHCSampler
    {
        class LHCPoint
        {
            public int x { get; set; }
            public int y { get; set; }
        }

        class OutPoint
        {
            public double x { get; set; }
            public double y { get; set; }
        }

        class TempSamplePointsForDimension
        {
            public List<double> points { get; private set; }
        }

        public class FinalSamplePoints
        {
            public double[] points { get; set; }
        }

        private int dimensions;
        private int samples;
        private List<LHCPoint> intervals;

        public LHCSampler(int dim, int samplesT)
        {
            init(dim, samplesT);
        }
        void init(int dim, int samplesT)
        {
            dimensions = dim;
            samples = samplesT;
            intervals = new List<LHCPoint>();
            
            for (int i = 0; i < samples; i++)
            {
                intervals.Add(new LHCPoint() {x = 0, y = 1});
            }
        }

        private TempSamplePointsForDimension divider(int dim)
        {
            TempSamplePointsForDimension ret = new TempSamplePointsForDimension();

            double dx = (double) (intervals[dim].y - intervals[dim].x) / samples;

            for (int i = 0; i < samples; i++)
            {
                double rnd = RNG.nextdouble();
                bool ok = (rnd >= 0) && (rnd <= 1);
                while (!ok)
                {
                    rnd = RNG.nextdouble();
                    ok = (rnd >= 0) && (rnd <= 1);
                }
                ret.points.Add(intervals[dim].x + (rnd * i) * dx);
            }

            return ret;
        }

        private void setInterval(int lower, int upper, int dim = 00)
        {
            if (dim != 0)
            {
                intervals[dim].x = lower;
                intervals[dim].y = upper;
            }
            else
            {
                for (int i = 0; i < dimensions; i++)
                {
                    intervals[dim].x = lower;
                }
            }
        }

        public FinalSamplePoints[]  getPoints(bool useShuffle = false)
        {
            List<TempSamplePointsForDimension> sampleset = new List<TempSamplePointsForDimension>();
            for (int dim = 0; dim < dimensions; dim++)
            {
                sampleset.Add(divider(dim));
            }

            FinalSamplePoints[] pointsForUse = new FinalSamplePoints[dimensions];
            for (int i = 0; i < dimensions; i++)
            {
                pointsForUse[i] = new FinalSamplePoints {points = new double[samples]};
            }

            if (useShuffle)
            {
                int[] indices = new int[samples];
                for (int i = 0; i < samples; i++)
                {
                    indices[i] = i;
                }
                for (int dimension = 0; dimension < dimensions; dimension++)
                {
                    indices.Shuffle();
                    for (int sample = 0; sample < samples; sample++)
                    {
                        pointsForUse[dimension].points[(samples - 1) - sample] = sampleset[sample].points[dimension];
                    }
                }
            }
            else
            {
                for (int sample = samples - 1; sample >= 0; sample--)
                {
                    for (int dimension = 0; dimension < dimensions; dimension++)
                    {
                        // We're pulling from the samples for each dimension as we shuffle, so need to recount each time.
                        int ptCount = sampleset[sample].points.Count;
                        int index = RNG.nextint(0, ptCount);
                        double val = sampleset[index].points[dimension];

                        // Remove our sample for this dimension
                        sampleset[sample].points.RemoveAt(dimension);
                        pointsForUse[dimension].points[(samples - 1) - sample] = val;
                    }
                }
            }

            return pointsForUse;
        }

        public void testLHC()
        {
            init(2,10);
            setInterval(-1, 1, 0);
            setInterval(-2, 1, 1);
            getPoints(true);
        }
    }
}