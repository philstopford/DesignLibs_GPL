using Eto.Drawing;

namespace VeldridEto;

public partial class VeldridDriver
{
	
	private SizeF ScreenToWorld(SizeF pt)
	{
		PointF pt1 = ScreenToWorld(0, 0);
		PointF pt2 = ScreenToWorld(pt.Width, pt.Height);
		return new SizeF(pt2.X - pt1.X, pt2.Y - pt1.Y);
	}

	private PointF ScreenToWorld(float x, float y)
	{
		// Use logical size for coordinate conversions to ensure consistent
		// behavior across different DPI settings
		double oX_2 = (double)Surface!.Width / 2;
		double oX_3 = ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oX_4 = x;

		double oXC = oX_4 - oX_2 + oX_3;


		return new PointF((x - (float)Surface.Width / 2) * (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom()) + ovpSettings.getCameraX(),
			((float)Surface.Height / 2 - y) * (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom()) + ovpSettings.getCameraY());
	}
	
	private Point WorldToScreen(float x, float y)
	{
		// int oX = (int)((x - ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom())) + Surface.Width / 2);

		// Use logical size for coordinate conversions to ensure consistent
		// behavior across different DPI settings
		double oX_2 = (double)Surface!.Width / 2;
		double oX_3 = ovpSettings.getCameraX() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oX_4 = x;

		int oXC = (int)(oX_4 - oX_3 + oX_2);

		// int oY = (int)((y - ovpSettings.getCameraY() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom())) + Surface.Height / 2);

		double oY_2 = (double)Surface.Height / 2;
		double oY_3 = ovpSettings.getCameraY() / (ovpSettings.getZoomFactor() * ovpSettings.getBaseZoom());
		double oY_4 = y;

		int oYC = (int)(oY_4 - oY_3 + oY_2);

		return new Point(oXC, oYC);
	}

	private Size WorldToScreen(SizeF pt)
	{
		Point pt1 = WorldToScreen(0, 0);
		Point pt2 = WorldToScreen(pt.Width, pt.Height);
		return new Size(pt2.X - pt1.X, pt2.Y - pt1.Y);
	}
}
