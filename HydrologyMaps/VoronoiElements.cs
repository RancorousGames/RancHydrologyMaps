/*
 * Created by SharpDevelop.
 * User: Burhan
 * Date: 17/06/2014
 * Time: 09:29 م
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */

/*
	  Copyright 2011 James Humphreys. All rights reserved.
	
	Redistribution and use in source and binary forms, with or without modification, are
	permitted provided that the following conditions are met:
	
	   1. Redistributions of source code must retain the above copyright notice, this list of
	      conditions and the following disclaimer.
	
	   2. Redistributions in binary form must reproduce the above copyright notice, this list
	      of conditions and the following disclaimer in the documentation and/or other materials
	      provided with the distribution.
	
	THIS SOFTWARE IS PROVIDED BY James Humphreys ``AS IS'' AND ANY EXPRESS OR IMPLIED
	WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
	FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
	ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	
	The views and conclusions contained in the software and documentation are those of the
	authors and should not be interpreted as representing official policies, either expressed
	or implied, of James Humphreys.
 */

/*
 * C# Version by Burhan Joukhadar
 * 
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */


using System;
using System.Collections.Generic;
using System.Drawing;
using HydrologyMaps;

namespace Voronoi2
{
	public class Point
	{
		public double X { get; set; }
		public double Y { get; set; }
	
		public Point ()
		{
		}

		public Point(double x, double y)
		{
			X = x;
			Y = y;
		}
		
		public double Distance(Point other)
		{
			double deltaX = other.X - this.X;
			double deltaY = other.Y - this.Y;
			return Math.Sqrt(deltaX * deltaX + deltaY * deltaY);
		}

		public static Point FromOtherPoint(HydrologyMaps.Point point)
		{
			return new Point(point.X, point.Y);
		}
		
		public static Point operator -(Point p1, Point p2)
		{
			return new Point(p1.X - p2.X, p1.Y - p2.Y);
		}
	}
	
	// use for sites and vertexes
	public class Site
	{
		public Point Coord { get; set; }
		public int SiteNumber { get; set; }

		public Site()
		{
			Coord = new Point(0, 0);
		}

		public Site(double x, double y)
		{
			Coord = new Point(x, y);
		}
	}
	
	public class Edge
	{
		public Edge(Site site1, Site site2)
		{
			EndPoint1 = site1;
			EndPoint2 = site2;
		}
		
		public Edge()
		{
		}


		public double A { get; set; } = 0;
		public double B { get; set; } = 0;
		public double C { get; set; } = 0;
		public Site EndPoint1 { get; set; }
		public Site EndPoint2 { get; set; }
		public Site Region1 { get; set; }
		public Site Region2 { get; set; }
		
		public int EdgeNumber { get; set; }
		
		
		
		public LineSegment ClipToBounds(double pxmin, double pymin, double pxmax, double pymax)
		{
			LineSegment clipped = new LineSegment();
    
			if (Math.Abs(A) > 1e-10)
			{
				clipped.Y1 = pymin;
				clipped.X1 = (C - B * clipped.Y1) / A;
        
				clipped.Y2 = pymax;
				clipped.X2 = (C - B * clipped.Y2) / A;
			}
			else
			{
				clipped.X1 = pxmin;
				clipped.Y1 = (C - A * clipped.X1) / B;

				clipped.X2 = pxmax;
				clipped.Y2 = (C - A * clipped.X2) / B;
			}

			if (A < 0 || (A == 0 && B > 0))
			{
				if (clipped.X1 < pxmin)
				{
					clipped.X1 = pxmin;
					clipped.Y1 = (C - A * clipped.X1) / B;
				}
				if (clipped.X2 > pxmax)
				{
					clipped.X2 = pxmax;
					clipped.Y2 = (C - A * clipped.X2) / B;
				}
			}
			else
			{
				if (clipped.X2 < pxmin)
				{
					clipped.X2 = pxmin;
					clipped.Y2 = (C - A * clipped.X2) / B;
				}
				if (clipped.X1 > pxmax)
				{
					clipped.X1 = pxmax;
					clipped.Y1 = (C - A * clipped.X1) / B;
				}
			}

			if (clipped.X1 < pxmin || clipped.X1 > pxmax || clipped.Y1 < pymin || clipped.Y1 > pymax ||
			    clipped.X2 < pxmin || clipped.X2 > pxmax || clipped.Y2 < pymin || clipped.Y2 > pymax)
			{
				return null;
			}

			return clipped;
		}
		
		
		
	}
	
	
	
	public class Halfedge
	{
		public Halfedge Left { get; set; }
		public Halfedge Right { get; set; }
		public Edge? Edge { get; set; }
		public bool Deleted { get; set; }
		public int PM { get; set; }
		public Site Vertex { get; set; }
		public double YStar { get; set; }
		public Halfedge Next { get; set; }

		public Halfedge()
		{
			Next = null;
		}
	}

	public class GraphEdge : IEdge
	{
		public GraphEdge(double x1, double y1, double x2, double y2)
		{
			X1 = x1;
			Y1 = y1;
			X2 = x2;
			Y2 = y2;
		}

		public GraphEdge(double x1, double y1, double x2, double y2, int leftSiteSiteNumber, int rightSiteSiteNumber)
		{
			X1 = x1;
			Y1 = y1;
			X2 = x2;
			Y2 = y2;
			Site1 = leftSiteSiteNumber;
			Site2 = rightSiteSiteNumber;
		}


		public double X1 { get; set; }
		public double Y1 { get; set; }
		public double X2 { get; set; }
		public double Y2 { get; set; }
		public int Site1 { get; set; }
		public int Site2 { get; set; }
		
		public (double x, double y) LeftSite => X1 < X2 ? (X1, Y1) : (X2, Y2);

		public (double x, double y) RightSite => X1 > X2 ? (X1, Y1) : (X2, Y2);

		public bool IsPointToLeft(Point point)
		{
			// Calculate the directional vector from the first coordinate to the point
			var dir = point - new Point(X1, Y1);

			// Calculate the directional vector of the edge
			var edgeDir = new Point(X2, Y2) - new Point(X1, Y1);

			// Calculate the cross product of the two vectors
			var crossProduct = dir.X * edgeDir.Y - dir.Y * edgeDir.X;

			// If the cross product is positive, the point is to the left of the edge
			return crossProduct > 0;
		}
		
		public Point? Intersection(GraphEdge other)
		{
			// Calculate the directional vectors for each edge
			var dir1 = new Point(X2, Y2) - new Point(X1, Y1);
			var dir2 = new Point(other.X2, other.Y2) - new Point(other.X1, other.Y1);

			// Calculate the determinant of the 2x2 matrix formed by the two directional vectors
			var det = dir1.X * dir2.Y - dir1.Y * dir2.X;

			// If the determinant is zero, the edges are parallel and do not intersect
			if (det == 0)
			{
				return null;
			}

			// Calculate the parameter values for the two edges
			var t1 = ((other.X1 - X1) * dir2.Y - (other.Y1 - Y1) * dir2.X) / det;
			var t2 = ((X1 - other.X1) * dir1.Y - (Y1 - other.Y1) * dir1.X) / -det;

			// If either of the parameter values is outside the range [0, 1], the edges do not intersect
			if (t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1)
			{
				return null;
			}

			// The edges intersect at the point given by the parameter value t1 (or t2)
			var x = X1 + t1 * dir1.X;
			var y = Y1 + t1 * dir1.Y;

			return new Point(x, y);
		}
	}
	
	// للترتيب
	public class SiteSorterYX : IComparer<Site>
	{
		public int Compare(Site p1, Site p2)
		{
			Point s1 = p1.Coord;
			Point s2 = p2.Coord;
			if (s1.Y < s2.Y) return -1;
			if (s1.Y > s2.Y) return 1;
			if (s1.X < s2.X) return -1;
			if (s1.X > s2.X) return 1;
			return 0;
		}
	}
	
	public class LineSegment
	{
		public double X1 { get; set; }
		public double Y1 { get; set; }
		public double X2 { get; set; }
		public double Y2 { get; set; }
	}
}
