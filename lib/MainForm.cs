/*
NPlot - A plotting library for .NET

Main.cs
Copyright (C) 2003-2004
Matt Howlett, Paolo Pierini

Port to Gtk# Miguel de icaza 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   
2. Redistributions in binary form must reproduce the following text in 
   the documentation and / or other materials provided with the 
   distribution: 
   
   "This product includes software developed as part of 
   the NPlot charting library project available from: 
   http://www.nplot.com/" 

------------------------------------------------------------------------

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Collections;
using System.ComponentModel;
using Gtk;
using System.Data;
using NPlot;
using System.IO;
using System.Reflection;


namespace NPlotDemo
{

	/// <summary>
	/// The main demo window.
	/// </summary>
	public class PlotSurface2DDemo : Gtk.Window
	{

		static void Main() 
		{
			Application.Init ();
			Window mf = new PlotSurface2DDemo ();
			mf.Show ();
			Application.Run ();
		}
		
		/// <summary>
		/// used to keep track of the current demo plot being displayed.
		/// </summary>
		private int currentPlot = 0;

		/// <summary>
		/// delegate for plot demo functions.
		/// </summary>
		private delegate void PlotDemoDelegate();

		/// <summary>
		///  list of the plot demos, initialized in the form constructor.
		/// </summary>
		private PlotDemoDelegate [] PlotRoutines;

		// Note that a NPlot.Windows.PlotSurface2D class
		// is used here. This has exactly the same 
		// functionality as the NPlot.PlotSurface2D 
		// class, except that it is derived from Forms.UserControl
		// and automatically paints itself in a windows.forms
		// application. Windows.PlotSurface2D can also paint itself
		// to other arbitrary Drawing.Graphics drawing surfaces
		// using the Draw method. (see printing later).

		Gtk.Button quitButton;
		Gtk.Button nextPlotButton;
		Gtk.Button prevPlotButton;
		NPlot.Gtk.PlotSurface2D plotSurface;

		double[] PlotQEExampleValues;
		string[] PlotQEExampleTextValues;

		public void PlotCircular()
		{
			this.plotSurface.Clear();
			plotSurface.Add( new HorizontalLine( 0.0, Color.LightGray ) );
			plotSurface.Add( new VerticalLine( 0.0, Color.LightGray ) );

			const int N = 400;
			const double start = -Math.PI * 7.0;
			const double end = Math.PI * 7.0;

			double[] xs = new double[N];
			double[] ys = new double[N];

			for (int i=0; i<N; ++i)
			{
				double t = ((double)i*(end - start)/(double)N + start);
				xs[i] = 0.5 * (t - 2.0 * Math.Sin(t));
				ys[i] = 2.0 * (1.0 - 2.0 * Math.Cos(t));
			}

			LinePlot lp = new LinePlot( ys, xs );
			lp.Pen = new Pen( Color.DarkBlue, 2.0f );
			lp.Label = "Circular Line"; // no legend, but still useful for copy data to clipboard.
			plotSurface.Add( lp );

			plotSurface.XAxis1 = new PiAxis( plotSurface.XAxis1 );

			plotSurface.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

			plotSurface.Refresh();
		}

		#region PlotWavelet
		public void PlotWavelet()
		{	
			this.plotSurface.Clear();

			// Create a new line plot from array data via the ArrayAdapter class.
			LinePlot lp = new LinePlot();
			lp.DataSource = makeDaub(256);
			lp.Color = Color.Green;
			lp.Label = "Daubechies Wavelet"; // no legend, but still useful for copy data to clipboard.

			Grid myGrid = new Grid();
			myGrid.VerticalGridType = Grid.GridType.Fine;
			myGrid.HorizontalGridType = Grid.GridType.Coarse;
			this.plotSurface.Add(myGrid);

			// And add it to the plot surface
			this.plotSurface.Add( lp );
			this.plotSurface.Title = "Reversed / Upside down Daubechies Wavelet";

			// Ok, the above will produce a decent default plot, but we would like to change
			// some of the Y Axis details. First, we'd like lots of small ticks (10) between 
			// large tick values. Secondly, we'd like to draw a grid for the Y values. To do 
			// this, we create a new LinearAxis (we could also use Label, Log etc). Rather than
			// starting from scratch, we use the constructor that takes an existing axis and
			// clones it (values in the superclass Axis only are cloned). PlotSurface2D
			// automatically determines a suitable axis when we add plots to it (merging
			// current requirements with old requirements), and we use this as our starting
			// point. Because we didn't specify which Y Axis we are using when we added the 
			// above line plot (there is one on the left - YAxis1 and one on the right - YAxis2)
			// PlotSurface2D.Add assumed we were using YAxis1. So, we create a new axis based on
			// YAxis1, update the details we want, then set the YAxis1 to be our updated one.
			LinearAxis myAxis = new LinearAxis( this.plotSurface.YAxis1 );
			myAxis.NumberOfSmallTicks = 9;
			this.plotSurface.YAxis1 = myAxis;
	
			// We would also like to modify the way in which the X Axis is printed. This time,
			// we'll just modify the relevant PlotSurface2D Axis directly. 
			this.plotSurface.XAxis1.WorldMax = 100.0f;
		
			this.plotSurface.PlotBackColor = Color.Beige;
			this.plotSurface.XAxis1.Reversed = true;
			this.plotSurface.YAxis1.Reversed = true;
		
			// Force a re-draw of the control. 
			this.plotSurface.Refresh();
		}
	

		private  float[] makeDaub( int len )
		{
			float[] daub4_h = 
			{ 0.482962913145f, 0.836516303737f, 0.224143868042f, -0.129409522551f };

			float[] daub4_g = 
			{ -0.129409522551f, -0.224143868042f, 0.836516303737f, -0.482962913145f };

			float[] a = new float[len];
			a[8] = 1.0f;
			float[] t;

			int ns = 4;  // number smooth
			while ( ns < len/2 ) 
			{
				t = (float[])a.Clone();

				ns *= 2;

				for ( int i=0; i<(ns*2); ++i ) 
				{
					a[i] = 0.0f;
				}

				// wavelet contribution
				for ( int i=0; i<ns; ++i ) 
				{
					for ( int j=0; j<4; ++j ) 
					{
						a[(2*i+j)%(2*ns)] += daub4_g[j] * t[i+ns];
					}
				}
				// smooth contribution
				for ( int i=0; i<ns; ++i ) 
				{
					for ( int j=0; j<4; ++j ) 
					{
						a[(2*i+j)%(2*ns)] += daub4_h[j]*t[i];
					}
				}
			}
			return a;
		}
		#endregion
		#region PlotLogAxis
		public void PlotLogAxis()
		{
			plotSurface.Clear();

			// draw a fine grid. 
			Grid fineGrid = new Grid();
			fineGrid.VerticalGridType = Grid.GridType.Fine;
			fineGrid.HorizontalGridType = Grid.GridType.Fine;
			plotSurface.Add( fineGrid );

			const int npt = 101;
			float[] x = new float[npt];
			float[] y = new float[npt];
			float step = 0.1f;
			for (int i=0; i<npt; ++i)
			{
				x[i] = i*step - 5.0f;
				y[i] = (float)Math.Pow( 10.0, x[i] );
			}
			float xmin = x[0];
			float xmax = x[npt-1];
			float ymin = (float)Math.Pow( 10.0, xmin );
			float ymax = (float)Math.Pow( 10.0, xmax );

			LinePlot lp = new LinePlot();
			lp.OrdinateData = y;
			lp.AbscissaData = x;
			lp.Pen = new Pen( Color.Red );
			plotSurface.Add( lp );

			LogAxis loga = new LogAxis( plotSurface.YAxis1 );
			loga.WorldMin = ymin;
			loga.WorldMax = ymax;
			loga.AxisColor = Color.Red;
			loga.LabelColor = Color.Red;
			loga.TickTextColor = Color.Red;
			loga.LargeTickStep = 1.0f;
			loga.Label = "10^x";
			plotSurface.YAxis1 = loga;

			LinePlot lp1 = new LinePlot();
			lp1.OrdinateData = y;
			lp1.AbscissaData = x;
			lp1.Pen = new Pen( Color.Blue );
			plotSurface.Add( lp1, PlotSurface2D.XAxisPosition.Bottom, PlotSurface2D.YAxisPosition.Right );
			LinearAxis lin = new LinearAxis( plotSurface.YAxis2 );
			lin.WorldMin = ymin;
			lin.WorldMax = ymax;
			lin.AxisColor = Color.Blue;
			lin.LabelColor = Color.Blue;
			lin.TickTextColor = Color.Blue;
			lin.Label = "10^x";
			plotSurface.YAxis2 = lin;
 
			LinearAxis lx = (LinearAxis)plotSurface.XAxis1;
			lx.WorldMin = xmin;
			lx.WorldMax = xmax;
			lx.Label = "x";

			//((LogAxis)plotSurface.YAxis1).LargeTickStep = 2;

			plotSurface.Title = "Mixed Linear/Log Axes";

			//plotSurface.XAxis1.LabelOffset = 20.0f;

			plotSurface.Refresh();
		}
		#endregion
		#region PlotLogLog
		public void PlotLogLog()
		{
			// log log plot
			plotSurface.Clear();

			Grid mygrid = new Grid();
			mygrid.HorizontalGridType = Grid.GridType.Fine;
			mygrid.VerticalGridType = Grid.GridType.Fine;
			plotSurface.Add(mygrid);

			int npt = 101;
			float [] x = new float[npt];
			float [] y = new float[npt];

			float step=0.1f;

			// plot a power law on the log-log scale
			for (int i=0; i<npt; ++i)
			{
				x[i] = (i+1)*step;
				y[i] = x[i]*x[i];
			}
			float xmin = x[0];
			float xmax = x[npt-1];
			float ymin = y[0];
			float ymax = y[npt-1];

			LinePlot lp = new LinePlot();
			lp.OrdinateData = y;
			lp.AbscissaData = x; 
			lp.Pen = new Pen( Color.Red );
			plotSurface.Add( lp );
			// axes
			// x axis
			LogAxis logax = new LogAxis( plotSurface.XAxis1 );
			logax.WorldMin = xmin;
			logax.WorldMax = xmax;
			logax.AxisColor = Color.Red;
			logax.LabelColor = Color.Red;
			logax.TickTextColor = Color.Red;
			logax.LargeTickStep = 1.0f;
			logax.Label = "x";
			plotSurface.XAxis1 = logax;
			// y axis
			LogAxis logay = new LogAxis( plotSurface.YAxis1 );
			logay.WorldMin = ymin;
			logay.WorldMax = ymax;
			logay.AxisColor = Color.Red;
			logay.LabelColor = Color.Red;
			logay.TickTextColor = Color.Red;
			logay.LargeTickStep = 1.0f;
			logay.Label = "x^2";
			plotSurface.YAxis1 = logay;

			LinePlot lp1 = new LinePlot();
			lp1.OrdinateData = y;
			lp1.AbscissaData = x;
			lp1.Pen = new Pen( Color.Blue );
			plotSurface.Add( lp1, PlotSurface2D.XAxisPosition.Top, PlotSurface2D.YAxisPosition.Right );
			// axes
			// x axis (lin)
			LinearAxis linx = (LinearAxis) plotSurface.XAxis2;
			linx.WorldMin = xmin;
			linx.WorldMax = xmax;
			linx.AxisColor = Color.Blue;
			linx.LabelColor = Color.Blue;
			linx.TickTextColor = Color.Blue;
			linx.Label = "x";
			plotSurface.XAxis2 = linx;
			// y axis (lin)
			LinearAxis liny = (LinearAxis) plotSurface.YAxis2;
			liny.WorldMin = ymin;
			liny.WorldMax = ymax;
			liny.AxisColor = Color.Blue;
			liny.LabelColor = Color.Blue;
			liny.TickTextColor = Color.Blue;
			liny.Label = "x^2";
			plotSurface.YAxis2 = liny;

			plotSurface.Title = "x^2 plotted with log(red)/linear(blue) axes";

			plotSurface.Refresh();
		}
		#endregion
		#region PlotSincFunction
		private void PlotSincFunction() 
		{
			plotSurface.Clear(); // clear everything. reset fonts. remove plot components etc.

			System.Random r = new Random();
			double[] a = new double[100];
			double[] b = new double[100];
			double mult = 0.00001f;
			for( int i=0; i<100; ++i )  
			{
				a[i] = ((double)r.Next(1000)/5000.0f-0.1f)*mult;
				if (i == 50 ) { b[i] = 1.0f*mult; } 
				else
				{
					b[i] = (double)Math.Sin((((double)i-50.0f)/4.0f))/(((double)i-50.0f)/4.0f);
					b[i] *= mult;
				}
				a[i] += b[i];
			}
		
			Marker m = new Marker(Marker.MarkerType.Cross1,6,new Pen(Color.Blue,2.0F));
			PointPlot pp = new PointPlot( m );
			pp.OrdinateData = a;
			pp.AbscissaData = new StartStep( -500.0, 10.0 );
			pp.Label = "Random";
			plotSurface.Add(pp); 

			LinePlot lp = new LinePlot();
			lp.OrdinateData = b;
			lp.AbscissaData = new StartStep( -500.0, 10.0 );
			lp.Pen = new Pen( Color.Red, 2.0f );
			plotSurface.Add( lp );

			plotSurface.Title = "Sinc Function";
			plotSurface.YAxis1.Label = "Magnitude";
			plotSurface.XAxis1.Label = "Position";

			Legend legend = new Legend();
			legend.AttachTo( PlotSurface2D.XAxisPosition.Top, PlotSurface2D.YAxisPosition.Left );
			legend.VerticalEdgePlacement = Legend.Placement.Inside;
			legend.HorizontalEdgePlacement = Legend.Placement.Inside;

			plotSurface.Legend = new Legend();

			plotSurface.Refresh();
		}
		#endregion
		#region PlotGaussian
		public void PlotGaussian()
		{
			plotSurface.Clear();
	
			System.Random r = new Random();
			
			int len = 35;
			double[] a = new double[len];
			double[] b = new double[len];

			for (int i=0; i<len; ++i) 
			{
				a[i] = (double)Math.Exp(-(double)(i-len/2)*(double)(i-len/2)/50.0f);
				b[i] = a[i] + (r.Next(10)/50.0f)-0.05f;
				if (b[i] < 0.0f) 
				{
					b[i] = 0;
				}
			}

			HistogramPlot sp = new HistogramPlot();
			sp.DataSource = b;
			sp.Pen = Pens.DarkBlue;
			sp.Filled = true;
			sp.RectangleBrush = new RectangleBrushes.HorizontalCenterFade( Color.Lavender, Color.Gold );
			sp.BaseWidth = 0.5f;
			sp.Label = "Random Data";
			LinePlot lp = new LinePlot();
			lp.DataSource = a;
			lp.Pen = new Pen( Color.Blue, 3.0f );
			lp.Label = "Gaussian Function";
			plotSurface.Add( sp );
			plotSurface.Add( lp );
			plotSurface.Legend = new Legend();
			plotSurface.YAxis1.WorldMin = 0.0f;
			plotSurface.Title = "Histogram Plot";
			plotSurface.Refresh();
		}
		#endregion
		#region PlotABC
		public void PlotABC()
		{
			plotSurface.Clear();
			const int size = 200;
			float [] xs = new float [size];
			float [] ys = new float [size];
			for (int i=0; i<size; i++)
			{
				xs[i] = (float)Math.Sin((double)i/(double)(size-1)*2.0*Math.PI);
				ys[i] = (float)Math.Cos((double)i/(double)(size-1)*6.0*Math.PI);
			}

			LinePlot lp = new LinePlot();
			lp.OrdinateData = ys;
			lp.AbscissaData = xs;
			Pen linePen = new Pen( Color.Yellow, 5.0f );
			lp.Pen = linePen;
			plotSurface.Add(lp);
			plotSurface.Title = "AxisConstraint.EqualScaling in action...";

			// Image downloaded from http://squidfingers.com. Thanks!
			Assembly a = Assembly.GetExecutingAssembly();
			System.IO.Stream file =
				a.GetManifestResourceStream( "NPlotDemo.resources.pattern01.jpg" );
			System.Drawing.Image im = System.Drawing.Image.FromStream( file );
			plotSurface.PlotBackImage = new Bitmap( im );

			plotSurface.AddAxesConstraint( new AxesConstraint.AspectRatio( 1.0, PlotSurface2D.XAxisPosition.Top, PlotSurface2D.YAxisPosition.Left ) );
			plotSurface.XAxis1.WorldMin = plotSurface.YAxis1.WorldMin;
			plotSurface.XAxis1.WorldMax = plotSurface.YAxis1.WorldMax;
			plotSurface.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

			plotSurface.Refresh();
		}
		#endregion
		#region PlotLabelAxis
		public void PlotLabelAxis()
		{
			plotSurface.Clear();

			Grid mygrid = new Grid();
			mygrid.VerticalGridType = Grid.GridType.Coarse;
			Pen majorGridPen = new Pen( Color.LightGray );
			float[] pattern = { 1.0f, 2.0f };
			majorGridPen.DashPattern = pattern;
			mygrid.MajorGridPen = majorGridPen;
			plotSurface.Add( mygrid );

			float[] xs = {20.0f, 31.0f, 27.0f, 38.0f, 24.0f, 3.0f, 2.0f };
			float[] xs2 = {7.0f, 10.0f, 42.0f, 9.0f, 2.0f, 79.0f, 70.0f };
			float[] xs3 = {1.0f, 20.0f, 20.0f, 25.0f, 10.0f, 30.0f, 30.0f };

			HistogramPlot hp = new HistogramPlot();
			hp.DataSource = xs;
			hp.BaseWidth = 0.6f;
			hp.RectangleBrush =
				new RectangleBrushes.HorizontalCenterFade( Color.FromArgb(255,255,200), Color.White );
			hp.Filled = true;
			hp.Label = "Developer Work";
			HistogramPlot hp2 = new HistogramPlot();
			hp2.DataSource = xs2;
			hp2.Label = "Web Browsing";
			hp2.RectangleBrush = RectangleBrushes.Horizontal.FaintGreenFade;
			hp2.Filled = true;
			hp2.StackedTo( hp );
			HistogramPlot hp3 = new HistogramPlot();
			hp3.DataSource = xs3;
			hp3.Label = "P2P Downloads";
			hp3.RectangleBrush = RectangleBrushes.Vertical.FaintBlueFade;
			hp3.Filled = true;
			hp3.StackedTo( hp2 );
			plotSurface.Add( hp );
			plotSurface.Add( hp2 );
			plotSurface.Add( hp3 );
			plotSurface.Legend = new Legend();

			LabelAxis la = new LabelAxis( plotSurface.XAxis1 );
			la.AddLabel( "Monday", 0.0f );
			la.AddLabel( "Tuesday", 1.0f );
			la.AddLabel( "Wednesday", 2.0f );
			la.AddLabel( "Thursday", 3.0f );
			la.AddLabel( "Friday", 4.0f );
			la.AddLabel( "Saturday", 5.0f );
			la.AddLabel( "Sunday", 6.0f );
			la.Label = "Days";
			la.TickTextFont = new Font( "Courier New", 8 );

			plotSurface.XAxis1 = la;
			plotSurface.YAxis1.WorldMin = 0.0;
			plotSurface.YAxis1.Label = "MBytes";
			((LinearAxis)plotSurface.YAxis1).NumberOfSmallTicks = 1;

			plotSurface.Title = "Internet useage for user johnc 09/01/03 - 09/07/03";

			plotSurface.XAxis1.TicksLabelAngle = 30.0f;

			plotSurface.PlotBackBrush = RectangleBrushes.Vertical.FaintRedFade;
			plotSurface.Refresh();
		}
		#endregion
		#region PlotParticles
		public void PlotParticles()
		{
			plotSurface.Clear();

			Grid mygrid = new Grid();
			mygrid.HorizontalGridType = Grid.GridType.Fine;
			mygrid.VerticalGridType = Grid.GridType.Fine;
			plotSurface.Add( mygrid );

			// in this example we synthetize a particle distribution
			// in the x-x' phase space and plot it, with the rms Twiss
			// ellipse and desnity distribution
			const int Particle_Number = 500;
			float [] x = new float[Particle_Number];
			float [] y = new float[Particle_Number];
			// Twiss parameters for the beam ellipse
			// 5 mm mrad max emittance, 1 mm beta function
			float alpha, beta, gamma, emit;
			alpha = -2.0f;
			beta = 1.0f;
			gamma = (1.0f + alpha * alpha) / beta;
			emit = 4.0f;

			float da, xmax, xpmax;
			da = -alpha / gamma;
			xmax = (float)Math.Sqrt(emit / gamma);
			xpmax = (float)Math.Sqrt(emit * gamma);

			Random rand = new Random();

			// cheap randomizer on the unit circle
			for (int i = 0; i<Particle_Number; i++)
			{
				float r;
				do
				{
					x[i] = (float)(2.0f * rand.NextDouble() - 1.0f);
					y[i] = (float)(2.0f * rand.NextDouble() - 1.0f);
					r = (float)Math.Sqrt(x[i] * x[i] + y[i] * y[i]);
				} while (r > 1.0f);
			}

			// transform to the tilted twiss ellipse
			for (int i =0; i<Particle_Number; ++i)
			{
				y[i] *= xpmax;
				x[i] = x[i] * xmax + y[i] * da;
			}
			plotSurface.Title = "Beam Horizontal Phase Space and Twiss ellipse";

			PointPlot pp = new PointPlot();
			pp.OrdinateData = y;
			pp.AbscissaData = x;
			pp.Marker = new Marker(Marker.MarkerType.FilledCircle ,4, new Pen(Color.Blue));
			plotSurface.Add(pp, PlotSurface2D.XAxisPosition.Bottom, PlotSurface2D.YAxisPosition.Left);

			// set axes
			LinearAxis lx = (LinearAxis) plotSurface.XAxis1;
			lx.Label = "Position - x [mm]";
			lx.NumberOfSmallTicks = 2;
			LinearAxis ly = (LinearAxis) plotSurface.YAxis1;
			ly.Label = "Divergence - x' [mrad]";
			ly.NumberOfSmallTicks = 2;
			
			// Draws the rms Twiss ellipse computed from the random data
			float [] xeli=new float [40];
			float [] yeli=new float [40];

			float a_rms, b_rms, g_rms, e_rms;

			Twiss(x, y, out a_rms, out b_rms, out g_rms, out e_rms);
			TwissEllipse(a_rms, b_rms, g_rms, e_rms, ref xeli, ref yeli);

			LinePlot lp = new LinePlot();
			lp.OrdinateData = yeli;
			lp.AbscissaData = xeli;
			plotSurface.Add(lp, PlotSurface2D.XAxisPosition.Bottom, PlotSurface2D.YAxisPosition.Left);
			lp.Pen = new Pen( Color.Red, 2.0f );
			// Draws the ellipse containing 100% of the particles
			// for a uniform distribution in 2D the area is 4 times the rms
			float [] xeli2 = new float [40];
			float [] yeli2 = new float [40];
			TwissEllipse(a_rms, b_rms, g_rms, 4.0F * e_rms, ref xeli2, ref yeli2);

			LinePlot lp2 = new LinePlot();
			lp2.OrdinateData = yeli2;
			lp2.AbscissaData = xeli2;
			plotSurface.Add( lp2, PlotSurface2D.XAxisPosition.Bottom, PlotSurface2D.YAxisPosition.Left );
			Pen p2 = new Pen( Color.Red, 2.0f );
			float [] pattern = { 5.0f, 40.0f };
			p2.DashPattern = pattern;
			lp2.Pen = p2;

			// now bin the particle position to create beam density histogram
			float range, min, max;
			min = (float)lx.WorldMin;
			max = (float)lx.WorldMax;
			range = max - min;

			const int Nbin = 30;
			float [] xbin = new float[Nbin+1];
			float [] xh = new float[Nbin+1];

			for (int j=0; j<=Nbin; ++j)
			{
				xbin[j] = min + j * range;
				if (j < Nbin) xh[j] = 0.0F;
			}
			for (int i =0; i<Particle_Number; ++i)
			{
				if (x[i] >= min && x[i] <= max)
				{
					int j;
					j = Convert.ToInt32(Nbin * (x[i] - min) / range);
					xh[j] += 1;
				}
			}
			StepPlot sp= new StepPlot();
			sp.OrdinateData = xh;
			sp.AbscissaData = new StartStep( min, range / Nbin );
			sp.Center = true;
			plotSurface.Add(sp, PlotSurface2D.XAxisPosition.Bottom, PlotSurface2D.YAxisPosition.Right);
			// axis formatting
			LinearAxis ly2 = (LinearAxis)plotSurface.YAxis2;
			ly2.WorldMin = 0.0f;
			ly2.Label = "Beam Density [a.u.]";
			ly2.NumberOfSmallTicks = 2;
			sp.Pen = new Pen( Color.Green, 2 );

			// Finally, refreshes the plot
			plotSurface.Refresh();
		}

		// Fill the array containing the rms twiss ellipse data points
		// ellipse is g*x^2+a*x*y+b*y^2=e
		private void TwissEllipse(float a, float b, float g, float e, ref float [] x, ref float [] y)
		{
			float rot, sr, cr, brot;
			if (a==0) 
			{
				rot=0;
			}
			else
			{
				rot=(float)(.5*Math.Atan(2.0 * a / (g - b)));
			}
			sr = (float)Math.Sin(rot);
			cr = (float)Math.Cos(rot);
			brot = g * sr * sr - 2.0F * a * sr * cr + b * cr * cr;
			int npt=x.Length;
			float theta;
		
			for (int i=0; i<npt;++i)
			{
				float xr,yr;
				theta = i * 2.0F * (float)Math.PI / (npt-1);
				xr = (float)(Math.Sqrt(e * brot) * Math.Cos(theta));
				yr = (float)(Math.Sqrt(e / brot) * Math.Sin(theta));
				x[i] = xr * cr - yr * sr;
				y[i] = xr * sr + yr * cr;
			}
		}
		// Evaluates the rms Twiss parameters from the particle coordinates
		private void Twiss(float [] x, float [] y, out float a, out float b, out float g, out float e)
		{
			float xave, xsqave, yave, ysqave, xyave;
			float sigmaxsq, sigmaysq, sigmaxy;
			int Npoints= x.Length;
			xave = 0;
			yave = 0;
			for (int i=0; i<Npoints; ++i)
			{
				xave += x[i];
				yave += y[i];
			}
			xave /= Npoints;
			yave /= Npoints;
			xsqave = 0;
			ysqave = 0;
			xyave = 0;
			for (int i=0;i<Npoints;i++)
			{
				xsqave += x[i] * x[i];
				ysqave += y[i] * y[i];
				xyave += x[i] * y[i];
			}
			xsqave /= Npoints;
			ysqave /= Npoints;
			xyave /= Npoints;
			sigmaxsq = xsqave - xave * xave;
			sigmaysq = ysqave - yave * yave;
			sigmaxy = xyave - xave * yave;
			// Now evaluates rms Twiss parameters
			e = (float)Math.Sqrt(sigmaxsq * sigmaysq - sigmaxy * sigmaxy);
			a = -sigmaxy / e;
			b = sigmaxsq / e;
			g = (1.0F + a * a) / b;
		}
		#endregion
		#region PlotQE
		public void PlotQE()
		{
			plotSurface.Clear();
			
			int len = 24;
			string[] s = new string[len];
			PlotQEExampleValues = new double[len];
			PlotQEExampleTextValues = new string[len];

			Random r = new Random();

			for (int i=0; i<len;i++)
			{
				PlotQEExampleValues[i] = 8.0f + 12.0f * (double)r.Next(10000) / 10000.0f;
				if (PlotQEExampleValues[i] > 18.0f)
				{
					PlotQEExampleTextValues[i] = "KCsTe";
				}
				else
				{
					PlotQEExampleTextValues[i] = "";
				}
				s[i] = i.ToString("00") + ".1";
			}

			PointPlot pp = new PointPlot();
			pp.DataSource = PlotQEExampleValues;
			pp.Marker = new Marker( Marker.MarkerType.Square, 10 );
			pp.Marker.DropLine = true;
			pp.Marker.Pen = Pens.CornflowerBlue;
			pp.Marker.Filled = false;
			plotSurface.Add( pp );

			LabelPointPlot tp1 = new LabelPointPlot();
			tp1.DataSource = PlotQEExampleValues;
			tp1.TextData = PlotQEExampleTextValues;
			tp1.LabelTextPosition = LabelPointPlot.LabelPositions.Above;
			tp1.Marker = new Marker( Marker.MarkerType.None, 10 );
			plotSurface.Add( tp1 );

			LabelAxis la = new LabelAxis( plotSurface.XAxis1 );
			for (int i=0; i<len; ++i)
			{
				la.AddLabel( s[i], i );
			}
			FontFamily ff = new FontFamily( "Verdana" );
			la.TickTextFont = new Font( ff, 7 );
			plotSurface.XAxis1 = la;

			plotSurface.Title = "Cs2Te Photocathode QE evolution";
			plotSurface.TitleFont = new Font(ff,15);
			plotSurface.XAxis1.WorldMin = -1.0f;
			plotSurface.XAxis1.WorldMax = len;
			plotSurface.XAxis1.LabelFont = new Font( ff, 10 );
			plotSurface.XAxis1.Label = "Cathode ID";
			plotSurface.YAxis1.Label = "QE [%]";
			plotSurface.YAxis1.LabelFont = new Font( ff, 10 );
			plotSurface.YAxis1.TickTextFont = new Font( ff, 10 );

			plotSurface.YAxis1.WorldMin = 0.0;
			plotSurface.YAxis1.WorldMax= 25.0;

			plotSurface.XAxis1.TicksLabelAngle = 60.0f;

			plotSurface.Refresh();
		}
		#endregion

		#region PlotDataSet
		void PlotDataSet()
		{
			plotSurface.Clear();

			// obtain stock information from xml file
			DataSet ds = new DataSet();
			System.IO.Stream file = 
				Assembly.GetExecutingAssembly().GetManifestResourceStream( "NPlotDemo.resources.asx_jbh.xml" );
			ds.ReadXml( file, System.Data.XmlReadMode.ReadSchema );
			DataTable dt = ds.Tables[0];

			// create CandlePlot.
			CandlePlot cp = new CandlePlot();
			cp.DataSource = dt;
			cp.AbscissaData = "Date";
			cp.OpenData = "Open";
			cp.LowData = "Low";
			cp.HighData = "High";
			cp.CloseData = "Close";
			cp.BearishColor = Color.Red;
			cp.BullishColor = Color.Green;
			cp.Style = CandlePlot.Styles.Filled;

			// calculate 10 day moving average and 2*sd line
			ArrayList av10 = new ArrayList();
			ArrayList sd2_10 = new ArrayList();
			ArrayList sd_2_10 = new ArrayList();
			ArrayList dates = new ArrayList();
			for (int i=0; i<dt.Rows.Count-10; ++i)
			{
				float sum = 0.0f;
				for (int j=0; j<10; ++j)
				{
					sum += (float)dt.Rows[i+j]["Close"];
				}
				float average = sum / 10.0f;
				av10.Add( average );
				sum = 0.0f;
				for (int j=0; j<10; ++j)
				{
					sum += ((float)dt.Rows[i+j]["Close"]-average)*((float)dt.Rows[i+j]["Close"]-average);
				}
				sum /= 10.0f;
				sum = 2.0f*(float)Math.Sqrt( sum );
				sd2_10.Add( average + sum );
				sd_2_10.Add( average - sum );
				dates.Add( (DateTime)dt.Rows[i+10]["Date"] );
			}

			// and a line plot of close values.
			LinePlot av = new LinePlot();
			av.OrdinateData = av10;
			av.AbscissaData = dates;
			av.Color = Color.LightGray;
			av.Pen.Width = 2.0f;

			LinePlot top = new LinePlot();
			top.OrdinateData = sd2_10;
			top.AbscissaData = dates;
			top.Color = Color.LightSteelBlue;
			top.Pen.Width = 2.0f;

			LinePlot bottom = new LinePlot();
			bottom.OrdinateData = sd_2_10;
			bottom.AbscissaData = dates;
			bottom.Color = Color.LightSteelBlue;
			bottom.Pen.Width = 2.0f;

			FilledRegion fr = new FilledRegion( top, bottom );
			//fr.RectangleBrush = new RectangleBrushes.Vertical( Color.FloralWhite, Color.GhostWhite );
			fr.RectangleBrush = new RectangleBrushes.Vertical( Color.FromArgb(255,255,240), Color.FromArgb(240,255,255) );
			plotSurface.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

			plotSurface.Add( fr );

			plotSurface.Add( new Grid() );

			plotSurface.Add( av );
			plotSurface.Add( top );
			plotSurface.Add( bottom );
			plotSurface.Add( cp );

			// now make an arrow... 
			ArrowItem arrow = new ArrowItem( new PointD( ((DateTime)dt.Rows[60]["Date"]).Ticks, 2.28 ), -80, "An interesting flat bit" );
			arrow.ArrowColor = Color.DarkBlue;
			arrow.PhysicalLength = 50;

			plotSurface.Add( arrow );

			plotSurface.Title = "AU:JBH";
			plotSurface.XAxis1.Label = "Date / Time";
			plotSurface.YAxis1.Label = "Price [$]";

			plotSurface.Refresh();

			plotSurface.Refresh();
		}
		#endregion
		#region PlotImage
		public void PlotImage()
		{
			// lazy way to simulate the readout of a text file containing a 19x19 map of datapoints
			string myfile = 
	"-1.251382E-3 -1.279191E-3 -7.230207E-4 -8.064462E-4 -5.005528E-4 -5.839783E-4 -1.696318E-3 -1.668509E-3 -3.893189E-4 -4.449358E-4 -1.473850E-3 -1.473850E-3 -1.974403E-3 -1.946594E-3 -2.085637E-3 -2.085637E-3 -1.612892E-3 -1.640701E-3 -1.863169E-3 " +
	"-1.251382E-3 -1.306999E-3 -6.674037E-4 -8.620631E-4 -4.449358E-4 -6.674037E-4 -1.668509E-4 -1.668509E-3 -3.615103E-4 -5.005528E-4 -5.561698E-5 -1.473850E-3 -4.449358E-4 -1.946594E-3 -6.395953E-4 -2.057828E-3 -1.585084E-3 -1.696318E-3 -1.807552E-3 " +
	"-1.223573E-3 -1.306999E-3 -6.117867E-4 -9.176802E-4 -3.893189E-4 -7.508292E-4 -1.390424E-4 -1.640701E-3 -3.615103E-4 -5.561698E-4 -5.561698E-5 -1.446041E-3 -4.449358E-4 -1.918786E-3 -6.117867E-4 -2.057828E-3 -1.585084E-3 -1.724126E-3 -1.779743E-3 " +
	"-1.251382E-3 -1.334807E-3 -5.839783E-4 -9.732971E-4 -3.615103E-4 -8.342547E-4 -1.390424E-4 3.893189E-4 1.751935E-3 2.919891E-3 3.476061E-3 3.031125E-3 1.807552E-3 6.674037E-4 -6.117867E-4 -2.030020E-3 -1.585084E-3 -1.779743E-3 -1.779743E-3 " +
    "-1.279191E-3 -1.362616E-3 -5.561698E-4 -1.028914E-3 -3.615103E-4 8.620631E-4 2.335913E-3 3.114551E-3 4.087848E-3 5.227996E-3 6.395952E-3 5.700740E-3 4.560592E-3 2.502764E-3 1.362616E-3 -6.117867E-4 -1.585084E-3 -1.807552E-3 -1.751935E-3 " +
    "-1.306999E-3 -1.390424E-3 -5.561698E-4 -1.056723E-3 1.890977E-3 4.087848E-3 6.117868E-3 9.621738E-3 1.357054E-2 1.721345E-2 1.715784E-2 1.462726E-2 1.059503E-2 6.368144E-3 3.253593E-3 1.279191E-3 6.674037E-4 -1.807552E-3 -1.779743E-3 " +
    "-1.390424E-3 -1.390424E-3 -5.561698E-4 1.585084E-3 4.560592E-3 8.481589E-3 1.437699E-2 2.155158E-2 2.702985E-2 3.078400E-2 3.134017E-2 2.892083E-2 2.338694E-2 1.446041E-2 6.757463E-3 3.031125E-3 6.674037E-4 -1.807552E-3 -1.807552E-3 " +
	"-1.446041E-3 -1.362616E-3 1.140148E-3 3.448253E-3 7.647335E-3 1.512782E-2 2.360941E-2 3.125674E-2 3.520555E-2 3.673501E-2 3.692967E-2 3.598418E-2 3.345361E-2 2.466613E-2 1.415452E-2 5.700740E-3 3.114551E-3 8.342547E-5 -1.835360E-3 " +
	"-1.529467E-3 -1.334807E-3 1.112340E-3 5.367038E-3 1.154052E-2 2.080075E-2 3.011659E-2 3.581733E-2 3.751365E-2 3.676282E-2 3.687406E-2 3.776393E-2 3.598418E-2 3.139579E-2 1.999430E-2 9.315844E-3 3.142359E-3 1.112340E-4 -1.863169E-3 " +
	"-1.640701E-3 -1.306999E-3 1.084531E-3 6.785271E-3 1.557275E-2 2.410996E-2 3.311991E-2 3.584514E-2 3.748584E-2 3.681844E-2 3.681844E-2 3.776393E-2 3.592857E-2 3.350923E-2 2.177405E-2 1.140148E-2 3.114551E-3 1.112340E-4 -1.890977E-3 " +
	"-1.696318E-3 -1.251382E-3 1.056723E-3 6.813080E-3 1.557275E-2 2.413777E-2 3.311991E-2 3.756927E-2 3.745804E-2 3.687406E-2 3.676282E-2 3.776393E-2 3.590076E-2 3.350923E-2 2.174624E-2 1.142929E-2 3.114551E-3 1.390424E-4 -1.918786E-3 " +
	"-1.779743E-3 -1.195765E-3 1.028914E-3 6.785271E-3 1.256944E-2 2.180186E-2 3.039468E-2 3.598418E-2 3.743023E-2 3.695748E-2 3.673501E-2 3.773612E-2 3.587295E-2 2.967166E-2 1.874292E-2 9.593928E-3 3.086742E-3 1.668509E-4 -1.918786E-3 " +
	"-1.863169E-3 -1.140148E-3 -1.362616E-3 3.058934E-3 7.285824E-3 1.532248E-2 2.472175E-2 3.195195E-2 3.478842E-2 3.701310E-2 3.670720E-2 3.545582E-2 3.220223E-2 2.333132E-2 1.354273E-2 5.978825E-3 1.279191E-3 1.668509E-4 -1.918786E-3 " +
	"-1.918786E-3 -1.084531E-3 -1.390424E-3 6.674037E-4 3.448253E-3 7.563909E-3 1.596207E-2 2.394311E-2 2.875398E-2 3.103427E-2 3.181291E-2 2.972727E-2 2.363721E-2 1.543371E-2 9.176801E-3 5.978825E-3 1.251382E-3 1.668509E-4 -1.890977E-3 " +
	"-1.974403E-3 -1.056723E-3 -1.390424E-3 6.674037E-4 1.362616E-3 3.781955E-3 8.036654E-3 1.304218E-2 1.824237E-2 2.127349E-2 2.174624E-2 1.963279E-2 1.573960E-2 1.148491E-2 7.341441E-3 3.197976E-3 -1.446041E-3 -1.557275E-3 -1.835360E-3 " +
	"-2.030020E-3 -1.056723E-3 -1.334807E-3 -1.084531E-3 -1.279191E-3 1.001106E-3 3.948805E-3 6.674037E-3 9.983247E-3 1.243039E-2 1.426576E-2 1.454384E-2 1.184642E-2 8.787482E-3 4.504975E-3 6.952122E-4 -1.418233E-3 -1.557275E-3 -1.779743E-3 " +
	"-2.057828E-3 -1.056723E-3 -1.306999E-3 -1.084531E-3 -1.251382E-3 -7.786377E-4 6.395953E-4 2.586189E-3 4.254699E-3 6.117868E-3 7.619526E-3 7.508292E-3 5.951017E-3 3.003317E-3 -9.176802E-4 -1.195765E-3 -1.362616E-3 -1.557275E-3 -1.751935E-3 " +
	"-2.085637E-3 -1.084531E-3 -1.251382E-3 -1.112340E-3 -1.223573E-3 -7.786377E-4 -8.620631E-4 -3.893189E-4 -4.449358E-4 3.337019E-4 1.167957E-3 9.454886E-4 3.893189E-4 -8.342547E-4 -9.176802E-4 -1.195765E-3 -1.306999E-3 -1.585084E-3 -1.696318E-3 " +
	"-2.113445E-3 -1.140148E-3 -1.195765E-3 -1.140148E-3 -1.167957E-3 -7.786377E-4 -8.342547E-4 -3.893189E-4 -4.171273E-4 -1.279191E-3 -1.251382E-3 -1.195765E-3 -1.195765E-3 -8.342547E-4 -8.620631E-4 -1.223573E-3 -1.251382E-3 -1.612892E-3 -1.668509E-3";
			string [] tokens = myfile.Split(new char [1] {' '});
			double [,] map = new double [19,19];
			for (int i=0; i < 19; ++i)
			{
				for (int j=0; j < 19; ++j)
				{
					map[i,j] = Convert.ToDouble(tokens[i*19+j], new
						System.Globalization.CultureInfo("en-US"));
				}
			}

			plotSurface.Clear();
			plotSurface.Title = "Cathode 11.2 QE Map";
			
			ImagePlot ip = new ImagePlot(map, -9.0f, 1.0f, -9.0f, 1.0f);

			plotSurface.Add(ip);
			plotSurface.XAxis1.WorldMin = -10.0f;
			plotSurface.XAxis1.WorldMax = 10.0f;
			plotSurface.XAxis1.Label = "x [mm]";
			//plotSurface.YAxis1.WorldMin = -10.0f;
			//plotSurface.YAxis1.WorldMax = 10.0f;
			plotSurface.YAxis1.Label = "y [mm]";

			plotSurface.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.None;

			//plotSurface.AddAxesConstraint( new AxesConstraint.AxisPosition( PlotSurface2D.YAxisPosition.Left, 0) );
			//plotSurface.AddAxesConstraint( new AxesConstraint.AxisPosition( PlotSurface2D.XAxisPosition.Top, 0.0f) );
			//plotSurface.AddAxesConstraint(
			//	new AxesConstraint.YPixelWorldLength(0.1f,PlotSurface2D.XAxisPosition.Bottom) );
			//plotSurface.AddAxesConstraint( new AxesConstraint.AspectRatio(1.0,PlotSurface2D.XAxisPosition.Top,PlotSurface2D.YAxisPosition.Left) );
			plotSurface.Refresh();

		}
		#endregion
		#region PlotMarkers
		public void PlotMarkers()
		{
		
			plotSurface.Clear();
			
			double[] y = new double[1] {1.0f};
			foreach (object i in Enum.GetValues(typeof(Marker.MarkerType)))
			{
				Marker m = new Marker( (Marker.MarkerType)Enum.Parse(typeof(Marker.MarkerType), i.ToString()), 8 );
				double[] x = new double[1];
				x[0] = (double) m.Type;
				PointPlot pp = new PointPlot();
				pp.OrdinateData = y;
				pp.AbscissaData = x;
				pp.Marker = m;
				pp.Label = m.Type.ToString();
				plotSurface.Add( pp );
			}
			plotSurface.Title = "Markers";
			plotSurface.YAxis1.Label = "Index";
			plotSurface.XAxis1.Label = "Marker";
			plotSurface.YAxis1.WorldMin = 0.0f;
			plotSurface.YAxis1.WorldMax = 2.0f;
			plotSurface.XAxis1.WorldMin -= 1.0f;
			plotSurface.XAxis1.WorldMax += 1.0f;

			Legend legend = new Legend();
			legend.AttachTo( PlotSurface2D.XAxisPosition.Top, PlotSurface2D.YAxisPosition.Right );
			legend.VerticalEdgePlacement = Legend.Placement.Outside;
			legend.HorizontalEdgePlacement = Legend.Placement.Inside;
			plotSurface.Legend = legend;

			plotSurface.Refresh();
		}
		#endregion
		#region PlotCandle
		void PlotCandle()
		{
			plotSurface.Clear();

			// obtain stock information from xml file
			DataSet ds = new DataSet();
			System.IO.Stream file =
				Assembly.GetExecutingAssembly().GetManifestResourceStream( "NPlotDemo.resources.asx_jbh.xml" );
			ds.ReadXml( file, System.Data.XmlReadMode.ReadSchema );
			DataTable dt = ds.Tables[0];

			// create CandlePlot.
			CandlePlot cp = new CandlePlot();
			cp.DataSource = dt;
			cp.AbscissaData = "Date";
			cp.OpenData = "Open";
			cp.LowData = "Low";
			cp.HighData = "High";
			cp.CloseData = "Close";
			cp.BearishColor = Color.Red;
			cp.BullishColor = Color.Green;
			cp.StickWidth = 3;
			cp.Color = Color.DarkBlue;

			plotSurface.Add( new Grid() );
			plotSurface.Add( cp );

			plotSurface.Title = "AU:JBH";
			plotSurface.XAxis1.Label = "Date / Time";
			plotSurface.YAxis1.Label = "Price [$]";
	
			plotSurface.Refresh();
		}
		#endregion
		#region PlotTest
		public void PlotTest()
		{

			plotSurface.Clear();

			// can plot different types.
			ArrayList l = new ArrayList();
			l.Add( (int)2 );
			l.Add( (double)1.0 );
			l.Add( (float)3.0 );
			l.Add( (int)5.0 );

			LinePlot lp1 = new LinePlot( new double[] {4.0, 3.0, 5.0, 8.0} );
			lp1.Pen = new Pen( Color.LightBlue );
			lp1.Pen.Width = 2.0f;

			//lp.AbscissaData = new StartStep( 0.0, 2.0 );

			LinePlot lp2 = new LinePlot( new double[] {2.0, 1.0, 4.0, 5.0} );
			lp2.Pen = new Pen( Color.LightBlue );
			lp2.Pen.Width = 2.0f;

			FilledRegion fr = new FilledRegion( lp1, lp2 );

			plotSurface.Add(fr);

			plotSurface.Add( new Grid() );
			plotSurface.Add(lp1);
			plotSurface.Add(lp2);

			ArrowItem a = new ArrowItem( new PointD( 2, 4 ), -50.0f, "Arrow" );
			a.HeadOffset = 5;
			a.ArrowColor = Color.Red;
			a.TextColor = Color.Purple;
			plotSurface.Add( a );

			MarkerItem m = new MarkerItem( new Marker( Marker.MarkerType.TriangleDown, 8, Color.ForestGreen ), 1.38, 2.9 );
			plotSurface.Add( m );

			plotSurface.XAxis1.TicksCrossAxis = true;
			
			((LinearAxis)plotSurface.XAxis1).LargeTickValue = -4.1;
			((LinearAxis)plotSurface.XAxis1).AutoScaleText = true;
			((LinearAxis)plotSurface.XAxis1).TicksIndependentOfPhysicalExtent = true;
			//plotSurface.XAxis1.Label = "Hello world";

			plotSurface.Refresh();

			/*
			plotSurface.AutoScaleTitle = false;
			plotSurface.AutoScaleAutoGeneratedAxes = true;
			
			plotSurface.Title = "My Title";

			double[] a = { 0, 2, 1, 4, double.NaN, double.NaN, 5, 8, 7, 9 };
			LinePlot lp = new LinePlot();
			lp.DataSource = a;
			lp.Label = "My Label";
			
			plotSurface.Add( lp );

			plotSurface.Legend = new Legend();
			plotSurface.Legend.AutoScaleText = false;
			plotSurface.Legend.NeverShiftAxes = true;
			plotSurface.Legend.HorizontalEdgePlacement = Legend.Placement.Inside;
			plotSurface.Legend.VerticalEdgePlacement = Legend.Placement.Inside;
			plotSurface.Legend.XOffset = -10;
			plotSurface.Legend.YOffset = 10;
			//plotSurface.AddAxesConstraint( new AxesConstraint.EqualSpacing() );

			((LinearAxis)plotSurface.XAxis1).Offset = 10.0;
			((LinearAxis)plotSurface.XAxis1).Scale = 27.0;
			//((LinearAxis)plotSurface.XAxis1).TicksIndependentOfPhysicalExtent = true;
			//((LinearAxis)plotSurface.YAxis1).TicksIndependentOfPhysicalExtent = true;

			AxesConstraint.AxisPosition c1 = 
				new NPlot.AxesConstraint.AxisPosition( PlotSurface2D.YAxisPosition.Left, 100.0f );

			AxesConstraint.AspectRatio c2 = 
				new AxesConstraint.AspectRatio( 5.0f, PlotSurface2D.YAxisPosition.Left );

			plotSurface.AddAxesConstraint( c1 );
			plotSurface.AddAxesConstraint( c2 );

			plotSurface.Refresh();
			*/
		}
		#endregion
		#region PlotWave
		public void PlotWave()
		{
			FileStream fs = new FileStream( @"c:\bounce.wav", System.IO.FileMode.Open );

			System.Int16[] w = new short[5000];
			byte[] a = new byte[10000];
			fs.Read( a, 0, 10000 );
			for (int i=0; i<5000; ++i)
			{
				w[i] = BitConverter.ToInt16(a,i*2);
			}

			plotSurface.Clear();

			LinePlot lp = new LinePlot();
			lp.DataSource = w;
			plotSurface.Add( lp );

			plotSurface.YAxis1.FlipTicksLabel = true;
			plotSurface.YAxis1.TicksLabelAngle = 200.0f;
			plotSurface.Refresh();

		}
		#endregion
		 
		public PlotSurface2DDemo() : base ("Gtk# NPlot Demo")
		{
			
			// List here the plot routines that you want to be accessed
			PlotRoutines = new PlotDemoDelegate [] {	//new PlotDemoDelegate(PlotWave),
														//new PlotDemoDelegate(PlotCandle),
				new PlotDemoDelegate(PlotCircular),
				//new PlotDemoDelegate(PlotDataSet),
													    new PlotDemoDelegate(PlotImage),
													    new PlotDemoDelegate(PlotQE),
													    new PlotDemoDelegate(PlotMarkers),
														new PlotDemoDelegate(PlotLogAxis),
														new PlotDemoDelegate(PlotLogLog),
													    new PlotDemoDelegate(PlotParticles), 
													    new PlotDemoDelegate(PlotWavelet), 
														new PlotDemoDelegate(PlotSincFunction), 
														new PlotDemoDelegate(PlotGaussian),
														new PlotDemoDelegate(PlotLabelAxis),
				//new PlotDemoDelegate(PlotABC)

														//new PlotDemoDelegate(PlotTest)
												};

			Setup ();
			
			// setup resize handler that takes care of placement of buttons, and sizing of
			// plotsurface2D when window is resized.

			// draw the first plot.
			currentPlot = 0;
			PlotRoutines[currentPlot]();
		}

		void Setup ()
		{
			quitButton = new Gtk.Button("Quit");
			quitButton.Clicked += quitButton_Click;
			
			nextPlotButton = new Gtk.Button("Next");
			nextPlotButton.Clicked += nextPlotButton_Click;
			
			prevPlotButton = new Gtk.Button("Prev");
			prevPlotButton.Clicked += prevPlotButton_Click;
			
			plotSurface = new NPlot.Gtk.PlotSurface2D();

			Gtk.Table t = new Gtk.Table (3, 3, false);
			Add (t);

			t.Attach (plotSurface, 0, 3, 0, 1);
			t.Attach (prevPlotButton, 0, 1, 1, 2, 0, 0, 0, 0);
			t.Attach (nextPlotButton, 1, 2, 1, 2, 0, 0, 0, 0);
			t.Attach (quitButton,     2, 3, 1, 2, 0, 0, 0, 0);
			ShowAll ();
		}

		/// <summary>
		/// callback for quit button click
		/// </summary>
		/// <param name="sender">unused</param>
		/// <param name="e">unused</param>
		private void quitButton_Click(object sender, System.EventArgs e)
		{
			Application.Quit ();
		}


		/// <summary>
		/// callback for next button click
		/// </summary>
		/// <param name="sender">unused</param>
		/// <param name="e">unused</param>
		private void nextPlotButton_Click(object sender, System.EventArgs e)
		{
			currentPlot += 1;
			if (currentPlot == PlotRoutines.Length)
			{
				currentPlot = 0;
			}

			int id = currentPlot+1;
			PlotRoutines[currentPlot]();
		}

		/// <summary>
		/// Callback for prev button click.
		/// </summary>
		/// <param name="sender">unused</param>
		/// <param name="e">unused</param>
		private void prevPlotButton_Click(object sender, System.EventArgs e)
		{
			currentPlot--;
			if( currentPlot == -1 ) currentPlot = PlotRoutines.Length-1;
			PlotRoutines[currentPlot]();
		}
	}
}
