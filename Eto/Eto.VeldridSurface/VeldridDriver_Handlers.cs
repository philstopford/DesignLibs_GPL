using Eto;
using Eto.Drawing;
using Eto.Forms;
using KDTree;

namespace VeldridEto;

public partial class VeldridDriver
{
	private void addKeyHandler(object sender, EventArgs e)
	{
		if (keyHandlerApplied)
		{
			return;
		}

		Surface.KeyDown += keyHandler;
		keyHandlerApplied = true;
	}

	private void removeKeyHandler(object sender, EventArgs e)
	{
		if (!keyHandlerApplied)
		{
			return;
		}

		hasFocus = false;
		Surface.KeyDown -= keyHandler;
		keyHandlerApplied = false;
	}

	private void keyHandler(object sender, KeyEventArgs e)
	{
		//e.Handled = true;

		if (ovpSettings.isLocked())
		{
			if (e.Key != Keys.F)
			{
				return;
			}
		}

		switch (e.Key)
		{
			case Keys.F:
				ovpSettings.lockVP(!ovpSettings.isLocked());
				break;
			case Keys.R:
				reset();
				break;
		}

		float stepping = 10.0f * ovpSettings.getZoomFactor();

		bool doUpdate = true;
		switch (e.Key)
		{
			case Keys.A:
				panHorizontal(-stepping);
				break;
			case Keys.D:
				panHorizontal(stepping);
				break;
			case Keys.W:
				panVertical(stepping);
				break;
			case Keys.S:
				panVertical(-stepping);
				break;
			case Keys.N:
				zoomOut(-1);
				break;
			case Keys.M:
				zoomIn(-1);
				break;
			case Keys.X:
				zoomExtents(-1);
				doUpdate = false; // update performed in extents
				break;
			case Keys.Z:
				zoomExtents(ovpSettings.selectedIndex);
				doUpdate = false; // update performed in extents
				break;
		}

		switch (doUpdate)
		{
			case true when Platform.Instance.IsGtk:
				updateHostFunc?.Invoke();
				break;
			case true:
				updateViewport();
				break;
		}
	}

	private void zoomHandler(object sender, MouseEventArgs e)
	{
		if (!ovpSettings.isLocked())
		{
			float wheelZoom = e.Delta.Height; // SystemInformation.MouseWheelScrollLines;
			switch (wheelZoom)
			{
				case > 0:
					zoomIn(wheelZoom);
					break;
				case < 0:
					zoomOut(-wheelZoom);
					break;
			}

			updateViewport();
		}

		e.Handled = true;
	}

	private void dragHandler(object sender, MouseEventArgs e)
	{
		if (!ovpSettings.isLocked())
		{
			switch (e.Buttons)
			{
				case MouseButtons.Primary:
				{
					PointF scaledLocation = e.Location * Surface.ParentWindow.LogicalPixelSize;

					switch (dragging)
					{
						case false:
							setDown(scaledLocation.X, scaledLocation.Y);
							break;
					}

					object locking = new();
					lock (locking)
					{
						float new_X = ovpSettings.getCameraX() - (scaledLocation.X - x_orig) *
							ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
						float new_Y = ovpSettings.getCameraY() + (scaledLocation.Y - y_orig) *
							ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom();
						ovpSettings.setCameraPos(new_X, new_Y);
						x_orig = scaledLocation.X;
						y_orig = scaledLocation.Y;
					}

					break;
				}
			}

			updateViewport();
		}

		e.Handled = true;
	}

	private void upHandler(object sender, MouseEventArgs e)
	{
		switch (e.Buttons)
		{
			case MouseButtons.Alternate:
			{
				menu?.Show(Surface);

				break;
			}
		}

		if (ovpSettings.isLocked())
		{
			return;
		}

		dragging = e.Buttons switch
		{
			MouseButtons.Primary => false,
			_ => dragging
		};
		e.Handled = true;
	}

	private void downHandler(object sender, MouseEventArgs e)
	{
		switch (e.Buttons)
		{
			case MouseButtons.Primary:
				setDown(e.Location.X, e.Location.Y);
				break;
		}

		if (e.Buttons == MouseButtons.Middle || e.Modifiers == Keys.Control && e.Buttons == MouseButtons.Primary)
		{
			selectByClick(e.Location.X, e.Location.Y);
		}

		e.Handled = true;
	}

	private void setFocus(object sender, EventArgs e)
	{
		if (hasFocus)
		{
			return;
		}

		Surface.Focus();
		hasFocus = true;
	}

	private void Clock_Elapsed(object sender, EventArgs e)
	{
		pUpdateViewport();
		if (done_drawing)
		{
			updateHostFunc?.Invoke();
			Surface.Invalidate();
			ovpSettings.changed = false;
			drawing = false;
			done_drawing = false;
		}
	}

	private void addKeyHandlers()
	{
		try
		{
			Surface.MouseDown += downHandler;
			Surface.MouseMove += dragHandler;
			Surface.MouseUp += upHandler;
			Surface.MouseWheel += zoomHandler;
			Surface.GotFocus += addKeyHandler;
			Surface.MouseEnter += addKeyHandler; // setFocus;
			Surface.LostFocus += removeKeyHandler;
			Surface.MouseLeave += removeKeyHandler;
		}
		catch (Exception ex)
		{
			Console.WriteLine("Ex: " + ex);
		}
	}

	private void setDown(float x, float y)
	{
		// might not be needed, but seemed like a safe approach to avoid re-setting these in a drag event.
		if (!dragging && !ovpSettings.isLocked())
		{
			x_orig = x;
			y_orig = y;
			dragging = true;
		}
	}

	private void selectByClick(float x, float y)
	{
		// Where did we click?
		PointF scaledLocation = new(x, y);
		scaledLocation = ScreenToWorld(scaledLocation.X * Surface.ParentWindow.LogicalPixelSize, scaledLocation.Y * Surface.ParentWindow.LogicalPixelSize);

		PointF cPos = ovpSettings.getCameraPos();

		// Populate our tree.
		int polyCount = ovpSettings.polyList.Count;
		switch (polyCount)
		{
			case > 0:
			{
				double[] distances = new double[polyCount];
				int[] indices = new int[polyCount];
				ParallelOptions po = new();
				//Parallel.For(0, polyCount, po, (poly, loopstate) =>
				for (int poly = 0; poly < ovpSettings.polyList.Count; poly++)
				{
					KDTree<PointF> pTree = new(2, ovpSettings.polyListPtCount[poly] + 1); // add one for the midpoint.
					foreach (PointF t1 in ovpSettings.polyList[poly].poly)
					{
						PointF t = new(t1.X, t1.Y);
						pTree.AddPoint(new double[] { t.X, t.Y }, t);
					}

					double maxX = ovpSettings.polyList[poly].poly.Max(p => p.X);
					double minX = ovpSettings.polyList[poly].poly.Min(p => p.X);
					double maxY = ovpSettings.polyList[poly].poly.Max(p => p.Y);
					double minY = ovpSettings.polyList[poly].poly.Min(p => p.Y);

					double deltaX = (maxX - minX) * 0.5f;
					double deltaY = (maxY - minY) * 0.5f;

					PointF midPoint = new((float)(minX + deltaX), (float)(minY + deltaY));
					pTree.AddPoint(new double[] { midPoint.X, midPoint.Y }, midPoint);

					// '1' forces a single nearest neighbor to be returned.
					NearestNeighbour<PointF> pIter = pTree.NearestNeighbors(new double[] { scaledLocation.X, scaledLocation.Y }, 1);
					while (pIter.MoveNext())
					{
						distances[poly] = Math.Abs(pIter.CurrentDistance);
						indices[poly] = ovpSettings.polySourceIndex[poly];
					}
				}
				//);

				int selIndex = indices[Array.IndexOf(distances, distances.Min())];

				updateHostSelectionFunc?.Invoke(selIndex);
				break;
			}
			default:
			{
				// Populate our tree.
				int lineCount = ovpSettings.lineList.Count;
				switch (lineCount)
				{
					case > 0:
					{
						double[] distances = new double[lineCount];
						int[] indices = new int[lineCount];
						ParallelOptions po = new();
						//Parallel.For(0, lineCount, po, (line, loopstate) =>
						for (int line = 0; line < ovpSettings.lineList.Count; line++)
						{
							KDTree<PointF> pTree = new(2, ovpSettings.lineListPtCount[line] + 1); // add one for the midpoint.
							foreach (PointF t1 in ovpSettings.lineList[line].poly)
							{
								PointF t = new(t1.X, t1.Y);
								pTree.AddPoint(new double[] { t.X, t.Y }, t);
							}

							double maxX = ovpSettings.lineList[line].poly.Max(p => p.X);
							double minX = ovpSettings.lineList[line].poly.Min(p => p.X);
							double maxY = ovpSettings.lineList[line].poly.Max(p => p.Y);
							double minY = ovpSettings.lineList[line].poly.Min(p => p.Y);

							double deltaX = (maxX - minX) * 0.5f;
							double deltaY = (maxY - minY) * 0.5f;

							PointF midPoint = new((float)(minX + deltaX), (float)(minY + deltaY));
							pTree.AddPoint(new double[] { midPoint.X, midPoint.Y }, midPoint);

							// '1' forces a single nearest neighbor to be returned.
							NearestNeighbour<PointF> pIter = pTree.NearestNeighbors(new double[] { scaledLocation.X, scaledLocation.Y }, 1);
							while (pIter.MoveNext())
							{
								distances[line] = Math.Abs(pIter.CurrentDistance);
								indices[line] = ovpSettings.lineSourceIndex[line];
							}
						}
						//);

						int selIndex = indices[Array.IndexOf(distances, distances.Min())];

						updateHostSelectionFunc?.Invoke(selIndex);
						break;
					}
				}

				break;
			}
		}
	}

	private float calcZoom(float delta)
	{
		float f = Math.Abs(delta) * 0.1f;
		f = delta switch
		{
			< 0 => 1.0f / f,
			_ => f
		};

		return f;
	}

	private void panVertical(float delta)
	{
		ovpSettings.setCameraY(ovpSettings.getCameraY() + delta / 10);
	}

	private void panHorizontal(float delta)
	{
		ovpSettings.setCameraX(ovpSettings.getCameraX() + delta / 10);
	}

	private void getExtents(int index)
	{
		if ((ovpSettings.polyList.Count == 0) && (ovpSettings.lineList.Count == 0))
		{
				ovpSettings.minX = 0;
				ovpSettings.maxX = 0;
				ovpSettings.minY = 0;
				ovpSettings.maxY = 0;
		}
		else
		{
			float minX = 0;
			float maxX = 0;
			float minY = 0, maxY = 0;

			bool set = false;

			if (ovpSettings.polyList.Count != 0)
			{
				for (int poly = 0; poly < ovpSettings.polyList.Count; poly++)
				{
					if (index != -1 && (ovpSettings.polySourceIndex[poly] != index || !ovpSettings.polyMask[poly]))
					{
						continue;
					}

					switch (set)
					{
						case false:
							minX = ovpSettings.polyList[poly].poly[0].X;
							maxX = ovpSettings.polyList[poly].poly[0].X;
							minY = ovpSettings.polyList[poly].poly[0].Y;
							maxY = ovpSettings.polyList[poly].poly[0].Y;
							set = true;
							break;
					}

					float tMinX = ovpSettings.polyList[poly].poly.Min(p => p.X);
					float tMaxX = ovpSettings.polyList[poly].poly.Max(p => p.X);
					float tMinY = ovpSettings.polyList[poly].poly.Min(p => p.Y);
					float tMaxY = ovpSettings.polyList[poly].poly.Max(p => p.Y);
					minX = Math.Min(minX, tMinX);
					maxX = Math.Max(maxX, tMaxX);
					minY = Math.Min(minY, tMinY);
					maxY = Math.Max(maxY, tMaxY);
				}
			}

			if (ovpSettings.lineList.Count != 0)
			{
				for (int line = 0; line < ovpSettings.lineList.Count; line++)
				{
					if (index != -1 && (ovpSettings.lineSourceIndex[line] != index || !ovpSettings.lineMask[line]))
					{
						continue;
					}

					switch (set)
					{
						case false:
							minX = ovpSettings.lineList[line].poly[0].X;
							maxX = ovpSettings.lineList[line].poly[0].X;
							minY = ovpSettings.lineList[line].poly[0].Y;
							maxY = ovpSettings.lineList[line].poly[0].Y;
							set = true;
							break;
					}

					float tMinX = ovpSettings.lineList[line].poly.Min(p => p.X);
					float tMaxX = ovpSettings.lineList[line].poly.Max(p => p.X);
					float tMinY = ovpSettings.lineList[line].poly.Min(p => p.Y);
					float tMaxY = ovpSettings.lineList[line].poly.Max(p => p.Y);
					minX = Math.Min(minX, tMinX);
					maxX = Math.Max(maxX, tMaxX);
					minY = Math.Min(minY, tMinY);
					maxY = Math.Max(maxY, tMaxY);
				}
			}
			ovpSettings.minX = minX;
			ovpSettings.maxX = maxX;
			ovpSettings.minY = minY;
			ovpSettings.maxY = maxY;
		}

		ovpSettings.changed = true;
	}
}
