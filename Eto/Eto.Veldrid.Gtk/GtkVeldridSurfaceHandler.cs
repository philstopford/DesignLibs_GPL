﻿using Eto.Drawing;
using Eto.GtkSharp.Forms;
using Eto.Veldrid;
using Eto.Veldrid.Gtk;
using Gtk;
using System;
using Veldrid;

[assembly: Eto.ExportHandler(typeof(VeldridSurface), typeof(GtkVeldridSurfaceHandler))]

namespace Eto.Veldrid.Gtk
{
	public class GtkVeldridSurfaceHandler : GtkControl<EtoEventBox, VeldridSurface, VeldridSurface.ICallback>, VeldridSurface.IHandler, VeldridSurface.IOpenGL
	{
		GLArea glArea;
		System.Action _makeCurrent;
		System.Action _clearCurrent;
		public Size RenderSize => Size.Round((SizeF)Widget.Size * Scale);

		float Scale => Widget.ParentWindow?.Screen?.LogicalPixelSize ?? 1;

		public override global::Gtk.Widget ContainerContentControl => glArea ?? base.ContainerContentControl;

		public GtkVeldridSurfaceHandler()
		{
			Control = new EtoEventBox { Handler = this };

			_makeCurrent = MakeCurrent;
			_clearCurrent = ClearCurrent;
		}

		public Swapchain CreateSwapchain()
		{
			Swapchain swapchain;

			if (Widget.Backend == GraphicsBackend.OpenGL)
			{
				swapchain = Widget.GraphicsDevice.MainSwapchain;
			}
			else
			{
				// To embed Veldrid in an Eto control, these platform-specific
				// versions of CreateSwapchain use the technique outlined here:
				//
				//   https://github.com/mellinoe/veldrid/issues/155
				//
				var source = SwapchainSource.CreateXlib(
					X11Interop.gdk_x11_display_get_xdisplay(Control.Display.Handle),
					X11Interop.gdk_x11_window_get_xid(Control.Window.Handle));

				var renderSize = RenderSize;
				swapchain = Widget.GraphicsDevice.ResourceFactory.CreateSwapchain(
					new SwapchainDescription(
						source,
						(uint)renderSize.Width,
						(uint)renderSize.Height,
						Widget.GraphicsDeviceOptions.SwapchainDepthFormat,
						Widget.GraphicsDeviceOptions.SyncToVerticalBlank,
						Widget.GraphicsDeviceOptions.SwapchainSrgbFormat));
			}

			return swapchain;
		}

		void glArea_InitializeGraphicsBackend(object sender, EventArgs e)
		{
			// Make context current to manually initialize a Veldrid GraphicsDevice.
			glArea.Context.MakeCurrent();
			Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
			// Veldrid clears the context at the end of initialization and sets it current in the worker thread.

			// Clear context in the worker thread for now to make Mesa happy.
			if (Widget.GraphicsDevice.GetOpenGLInfo(out BackendInfoOpenGL glInfo))
			{
				// This action has to wait so GTK can manage the context after this method.
				glInfo.ExecuteOnGLThread(_clearCurrent);//, wait: true);
			}

			glArea.Render += glArea_Render;
			glArea.Resize += glArea_Resize;
		}

		void Control_InitializeGraphicsBackend(object sender, EventArgs e)
		{
			Callback.OnInitializeBackend(Widget, new InitializeEventArgs(RenderSize));
		}

		bool skipDraw;

		private void glArea_Resize(object o, ResizeArgs args)
		{
			skipDraw = false;
			Callback.OnResize(Widget, new ResizeEventArgs(RenderSize));
		}

		void glArea_Render(object o, RenderArgs args)
		{
			if (!skipDraw)
			{
				skipDraw = true;

				// GTK makes the context current for us, so we need to clear it to hand it over to the Veldrid worker.
				Gdk.GLContext.ClearCurrent();

				// Make context current on the Veldrid worker.
				if (Widget.GraphicsDevice.GetOpenGLInfo(out BackendInfoOpenGL glInfo))
				{
					// No need for this action to wait, we just need it done some time before issuing commands.
					glInfo.ExecuteOnGLThread(_makeCurrent);//, wait: false);
				}

				// It's important to only issue Veldrid commands in OnDraw,
				// since we only have a GL context current in the worker here.
				Callback.OnDraw(Widget, EventArgs.Empty);

				// Clear the context from the worker so GTK can use it again.
				// This action needs to wait so the context is cleared before we continue
				// out of this method (either for the coming MakeCurrent or for GTK).
				glInfo?.ExecuteOnGLThread(_clearCurrent);//, wait: true);

				// GTK does not seem to need the context current after the Render event,
				// but setting it back is safer if this assumption changes in the future.
				glArea?.MakeCurrent();
			}
			skipDraw = false;
		}

		private void MakeCurrent()
		{
			glArea?.MakeCurrent();
		}

		private void ClearCurrent()
		{
			Gdk.GLContext.ClearCurrent();
		}

		// TODO: Figure this one out! The docstring for this property in Veldrid's OpenGLPlatformInfo is ambiguous.
		IntPtr VeldridSurface.IOpenGL.OpenGLContextHandle => glArea?.Context.Handle ?? IntPtr.Zero;

		IntPtr VeldridSurface.IOpenGL.GetProcAddress(string name) => X11Interop.glXGetProcAddress(name);

		void VeldridSurface.IOpenGL.MakeCurrent(IntPtr context) => MakeCurrent();

		IntPtr VeldridSurface.IOpenGL.GetCurrentContext() => Gdk.GLContext.Current.Handle;

		void VeldridSurface.IOpenGL.ClearCurrentContext() => ClearCurrent();

		void VeldridSurface.IOpenGL.DeleteContext(IntPtr context)
		{
		}

		void VeldridSurface.IOpenGL.SwapBuffers()
		{
			// GLArea doesn't support drawing directly, so we queue a render but don't actually call OnDraw
			if (skipDraw)
				return;

			skipDraw = true;
			glArea?.QueueRender();
		}

		void VeldridSurface.IOpenGL.SetSyncToVerticalBlank(bool on)
		{
		}

		void VeldridSurface.IOpenGL.SetSwapchainFramebuffer()
		{
		}

		void VeldridSurface.IOpenGL.ResizeSwapchain(uint width, uint height)
		{
		}

		void Eto.Forms.Control.IHandler.Invalidate(Rectangle rect, bool invalidateChildren)
		{
			skipDraw = false;
			glArea?.QueueRender();
		}

		void Eto.Forms.Control.IHandler.Invalidate(bool invalidateChildren)
		{
			skipDraw = false;
			glArea?.QueueRender();
		}


		protected override void Initialize()
		{
			base.Initialize();

			if (Widget.Backend == GraphicsBackend.OpenGL)
			{
				glArea = new GLArea();
				glArea.CanFocus = true;

				// Veldrid technically supports as low as OpenGL 3.0, but the full
				// complement of features is only available with 3.3 and higher.
				glArea.SetRequiredVersion(3, 3);

				glArea.HasDepthBuffer = true;
				glArea.HasStencilBuffer = true;
				Control.Child = glArea;
				glArea.Realized += glArea_InitializeGraphicsBackend;
			}
			else
			{
				Control.CanFocus = true;
				Control.Realized += Control_InitializeGraphicsBackend;
			}

		}
	}
}
