using System;
using System.Runtime.CompilerServices;

namespace KDTree;

/// <summary>
/// A KD-Tree node which supports a generic number of dimensions.  All data items
/// need the same number of dimensions.
/// This node splits based on the largest range of any dimension.
/// </summary>
/// <typeparam name="T">The generic data type this structure contains.</typeparam>
/// <remarks>This is based on this: https://bitbucket.org/rednaxela/knn-benchmark/src/tip/ags/utils/dataStructures/trees/thirdGenKD/ </remarks>
public class KDNode<T>
{
    #region Internal properties and constructor
    // All types
    /// <summary>
    /// The number of dimensions for this node.
    /// </summary>
    protected internal int iDimensions;

    /// <summary>
    /// The maximum capacity of this node.
    /// </summary>
    protected internal int iBucketCapacity;

    // Leaf only
    /// <summary>
    /// The array of locations.  [index][dimension]
    /// </summary>
    protected internal double[][] tPoints;

    /// <summary>
    /// The array of data values. [index]
    /// </summary>
    protected internal T[] tData;

    // Stem only
    /// <summary>
    /// The left and right children.
    /// </summary>
    protected internal KDNode<T> pLeft, pRight;
    /// <summary>
    /// The split dimension.
    /// </summary>
    protected internal int iSplitDimension;
    /// <summary>
    /// The split value (larger go into the right, smaller go into left)
    /// </summary>
    protected internal double fSplitValue;

    // Bounds
    /// <summary>
    /// The min and max bound for this node.  All dimensions.
    /// </summary>
    protected internal double[] tMinBound, tMaxBound;

    /// <summary>
    /// Does this node represent only one point.
    /// </summary>
    protected internal bool bSinglePoint;

    /// <summary>
    /// Protected method which constructs a new KDNode.
    /// </summary>
    /// <param name="iDimensions">The number of dimensions for this node (all the same in the tree).</param>
    /// <param name="iBucketCapacity">The initial capacity of the bucket.</param>
    protected KDNode(int iDimensions, int iBucketCapacity)
    {
        // Variables.
        this.iDimensions = iDimensions;
        this.iBucketCapacity = iBucketCapacity;
        Size = 0;
        bSinglePoint = true;

        // Setup leaf elements.
        tPoints = new double[iBucketCapacity + 1][];
        tData = new T[iBucketCapacity + 1];
    }
    #endregion

    #region External Operations
    /// <summary>
    /// The number of items in this leaf node and all children.
    /// </summary>
    public int Size { get; private set; }

    /// <summary>
    /// Is this KDNode a leaf or not?
    /// </summary>
    public bool IsLeaf => tPoints != null;

    /// <summary>
    /// Insert a new point into this leaf node.
    /// </summary>
    /// <param name="tPoint">The position which represents the data.</param>
    /// <param name="kValue">The value of the data.</param>
    public void AddPoint(double[] tPoint, T kValue)
    {
        // Find the correct leaf node.
        KDNode<T> pCursor = this;
        while (!pCursor.IsLeaf)
        {
            // Extend the size of the leaf.
            pCursor.ExtendBounds(tPoint);
            pCursor.Size++;

            // If it is larger select the right, or lower,  select the left.
            pCursor = tPoint[pCursor.iSplitDimension] > pCursor.fSplitValue ? pCursor.pRight : pCursor.pLeft;
        }

        // Insert it into the leaf.
        pCursor.AddLeafPoint(tPoint, kValue);
    }
    #endregion

    #region Internal Operations
    /// <summary>
    /// Insert the point into the leaf.
    /// </summary>
    /// <param name="tPoint">The point to insert the data at.</param>
    /// <param name="kValue">The value at the point.</param>
    private void AddLeafPoint(double[] tPoint, T kValue)
    {
        // Add the data point to this node.
        tPoints[Size] = tPoint;
        tData[Size] = kValue;
        ExtendBounds(tPoint);
        Size++;

        // Split if the node is getting too large in terms of data.
        if (Size != tPoints.Length - 1)
        {
            return;
        }

        // If the node is getting too physically large.
        if (CalculateSplit())
        {
            // If the node successfully had it's split value calculated, split node.
            SplitLeafNode();
        }
        else
        {
            // If the node could not be split, enlarge node data capacity.
            IncreaseLeafCapacity();
        }
    }

    /// <summary>
    /// If the point lies outside the boundaries, return false else true.
    /// </summary>
    /// <param name="tPoint">The point.</param>
    /// <returns>True if the point is inside the boundaries, false outside.</returns>
    private bool CheckBounds(double[] tPoint)
    {
        for (int i = 0; i < iDimensions; ++i)
        {
            if (tPoint[i] > tMaxBound[i])
            {
                return false;
            }

            if (tPoint[i] < tMinBound[i])
            {
                return false;
            }
        }
        return true;
    }

    /// <summary>
    /// Extend this node to contain a new point.
    /// </summary>
    /// <param name="tPoint">The point to contain.</param>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void ExtendBounds(double[] tPoint)
    {
        // If we don't have bounds, create them using the new point then bail.
        if (tMinBound == null)
        {
            tMinBound = new double[iDimensions];
            tMaxBound = new double[iDimensions];
            // Use Span for better performance than Array.Copy
            tPoint.AsSpan(0, iDimensions).CopyTo(tMinBound.AsSpan());
            tPoint.AsSpan(0, iDimensions).CopyTo(tMaxBound.AsSpan());
            return;
        }

        // For each dimension.
        for (int i = 0; i < iDimensions; ++i)
        {
            switch (tPoint[i])
            {
                case double.NaN:
                {
                    if (!double.IsNaN(tMinBound[i]) || !double.IsNaN(tMaxBound[i]))
                    {
                        bSinglePoint = false;
                    }

                    tMinBound[i] = double.NaN;
                    tMaxBound[i] = double.NaN;
                    break;
                }
                default:
                {
                    if (tMinBound[i] > tPoint[i])
                    {
                        tMinBound[i] = tPoint[i];
                        bSinglePoint = false;
                    }
                    else if (tMaxBound[i] < tPoint[i])
                    {
                        tMaxBound[i] = tPoint[i];
                        bSinglePoint = false;
                    }

                    break;
                }
            }
        }
    }

    /// <summary>
    /// Double the capacity of this leaf.
    /// </summary>
    private void IncreaseLeafCapacity()
    {
        Array.Resize<double[]>(ref tPoints, tPoints.Length * 2);
        Array.Resize<T>(ref tData, tData.Length * 2);
    }

    /// <summary>
    /// Work out if this leaf node should split.  If it should, a new split value and dimension is calculated
    /// based on the dimension with the largest range.
    /// </summary>
    /// <returns>True if the node split, false if not.</returns>
    private bool CalculateSplit()
    {
        switch (bSinglePoint)
        {
            // Don't split if we are just one point.
            case true:
                return false;
        }

        // Find the dimension with the largest range.  This will be our split dimension.
        double fWidth = 0;
        for (int i = 0; i < iDimensions; i++)
        {
            double fDelta = tMaxBound[i] - tMinBound[i];
            fDelta = fDelta switch
            {
                double.NaN => 0,
                _ => fDelta
            };

            if (!(fDelta > fWidth))
            {
                continue;
            }

            iSplitDimension = i;
            fWidth = fDelta;
        }

        switch (fWidth)
        {
            // If we are not wide (i.e. all the points are in one place), don't split.
            case 0:
                return false;
        }

        fSplitValue = fSplitValue switch
        {
            // Never split on infinity or NaN.
            double.PositiveInfinity => double.MaxValue,
            double.NegativeInfinity => double.MinValue,
            // Split in the middle of the node along the widest dimension.
            _ => (tMinBound[iSplitDimension] + tMaxBound[iSplitDimension]) * 0.5
        };

        // Don't let the split value be the same as the upper value as
        // can happen due to rounding errors!
        if (Math.Abs(fSplitValue - tMaxBound[iSplitDimension]) <= double.Epsilon)
        {
            fSplitValue = tMinBound[iSplitDimension];
        }

        // Success
        return true;
    }

    /// <summary>
    /// Split this leaf node by creating left and right children, then moving all the children of
    /// this node into the respective buckets.
    /// </summary>
    private void SplitLeafNode()
    {
        // Create the new children.
        pRight = new KDNode<T>(iDimensions, iBucketCapacity);
        pLeft = new KDNode<T>(iDimensions, iBucketCapacity);

        // Move each item in this leaf into the children.
        for (int i = 0; i < Size; ++i)
        {
            // Store.
            double[] tOldPoint = tPoints[i];
            T kOldData = tData[i];

            // If larger, put it in the right.
            if (tOldPoint[iSplitDimension] > fSplitValue)
            {
                pRight.AddLeafPoint(tOldPoint, kOldData);
            }

            // If smaller, put it in the left.
            else
            {
                pLeft.AddLeafPoint(tOldPoint, kOldData);
            }
        }

        // Wipe the data from this KDNode.
        tPoints = null;
        tData = null;
    }
    #endregion
}