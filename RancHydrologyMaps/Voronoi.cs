/*
 * Created by SharpDevelop.
 * User: Burhan
 * Date: 17/06/2014
 * Time: 11:30 م
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */

/*
* The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
* Bell Laboratories.
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

/* 
 * This code was originally written by Stephan Fortune in C code.  I, Shane O'Sullivan,
 * have since modified it, encapsulating it in a C++ class and, fixing memory leaks and
 * adding accessors to the Voronoi Edges.
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

/* 
 * Java Version by Zhenyu Pan
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
using System.Net;
using HydrologyMaps;
using NetTopologySuite.Geometries;

namespace Voronoi2
{
    /// <summary>
    /// Description of Voronoi.
    /// </summary>
    public class Voronoi
    {
        // ************* Private members ******************
        double borderMinX, borderMaxX, borderMinY, borderMaxY;
        int siteidx;
        double xmin, xmax, ymin, ymax, deltax, deltay;
        int nvertices;
        int nedges;
        int nsites;
        Site[] sites;
        Site bottomsite;
        int sqrt_nsites;
        double minDistanceBetweenSites;
        int PQcount;
        int PQmin;
        int PQhashsize;
        Halfedge[] PQhash;

        const int LE = 0;
        const int RE = 1;

        int ELhashsize;
        Halfedge[] ELhash;
        Halfedge ELleftend, ELrightend;
        List<GraphEdge> allEdges;

        // ************* Public methods ******************
        // ******************************************

        // constructor
        public Voronoi(double minDistanceBetweenSites)
        {
            siteidx = 0;
            sites = null;

            allEdges = null;
            this.minDistanceBetweenSites = minDistanceBetweenSites;
        }

        /**
		 * 
		 * @param minX The minimum X of the bounding box around the voronoi
		 * @param maxX The maximum X of the bounding box around the voronoi
		 * @param minY The minimum Y of the bounding box around the voronoi
		 * @param maxY The maximum Y of the bounding box around the voronoi
		 * @return
		 */
        public List<GraphEdge> GenerateVoronoi(List<IGraphNode> graphNodes, double minX, double maxX, double minY,
            double maxY)
        {
            Sort(graphNodes);

            // Check bounding box inputs - if mins are bigger than maxes, swap them
            double temp = 0;
            if (minX > maxX)
            {
                temp = minX;
                minX = maxX;
                maxX = temp;
            }

            if (minY > maxY)
            {
                temp = minY;
                minY = maxY;
                maxY = temp;
            }

            borderMinX = minX;
            borderMinY = minY;
            borderMaxX = maxX;
            borderMaxY = maxY;

            siteidx = 0;
            VoronoiBd();
            return allEdges;
        }


        /*********************************************************
         * Private methods - implementation details
         ********************************************************/

        private void Sort(List<IGraphNode> graphNodes)
        {
            sites = null;
            allEdges = new List<GraphEdge>();

            nsites = graphNodes.Count;
            nvertices = 0;
            nedges = 0;

            double sn = (double)nsites + 4;
            sqrt_nsites = (int)Math.Sqrt(sn);

            sites = new Site[nsites];
            xmin = graphNodes[0].Point.X;
            ymin = graphNodes[0].Point.Y;
            xmax = graphNodes[0].Point.X;
            ymax = graphNodes[0].Point.Y;

            for (int i = 0; i < nsites; i++)
            {
                sites[i] = new Site();
                sites[i].Coord = Point.FromOtherPoint(graphNodes[i].Point);
                sites[i].SiteNumber = i;

                if (graphNodes[i].Point.X < xmin)
                    xmin = graphNodes[i].Point.X;
                else if (graphNodes[i].Point.X > xmax)
                    xmax = graphNodes[i].Point.X;

                if (graphNodes[i].Point.Y < ymin)
                    ymin = graphNodes[i].Point.Y;
                else if (graphNodes[i].Point.Y > ymax)
                    ymax = graphNodes[i].Point.Y;
            }

            QSort(sites);
            deltax = xmax - xmin;
            deltay = ymax - ymin;
        }

        private void QSort(Site[] sitesToSort)
        {
            List<Site> siteList = new List<Site>(sitesToSort.Length);
            for (int i = 0; i < sitesToSort.Length; i++)
            {
                siteList.Add(sitesToSort[i]);
            }

            siteList.Sort(new SiteSorterYX());

            // Copy back into the array
            for (int i = 0; i < sitesToSort.Length; i++)
            {
                sitesToSort[i] = siteList[i];
            }
        }

        private void SortNode(List<IGraphNode> graphNodes)
        {
            nsites = graphNodes.Count;
            sites = new Site[nsites];
            xmin = graphNodes[0].Point.X;
            ymin = graphNodes[0].Point.Y;
            xmax = graphNodes[0].Point.X;
            ymax = graphNodes[0].Point.Y;

            for (int i = 0; i < nsites; i++)
            {
                sites[i] = new Site();
                sites[i].Coord = new Point(graphNodes[i].Point.X, graphNodes[i].Point.Y);
                sites[i].SiteNumber = i;

                if (graphNodes[i].Point.X < xmin)
                    xmin = graphNodes[i].Point.X;
                else if (graphNodes[i].Point.X > xmax)
                    xmax = graphNodes[i].Point.X;

                if (graphNodes[i].Point.Y < ymin)
                    ymin = graphNodes[i].Point.Y;
                else if (graphNodes[i].Point.Y > ymax)
                    ymax = graphNodes[i].Point.Y;
            }

            QSort(sites);
            deltax = xmax - xmin;
            deltay = ymax - ymin;
        }


        private Site NextOne()
        {
            Site s;
            if (siteidx < nsites)
            {
                s = sites[siteidx];
                siteidx++;
                return s;
            }

            return null;
        }

        private Edge Bisect(Site s1, Site s2)
        {
            double dx, dy, adx, ady;
            Edge newEdge;

            newEdge = new Edge();

            newEdge.Region1 = s1;
            newEdge.Region2 = s2;

            newEdge.EndPoint1 = null;
            newEdge.EndPoint2 = null;

            dx = s2.Coord.X - s1.Coord.X;
            dy = s2.Coord.Y - s1.Coord.Y;

            adx = dx > 0 ? dx : -dx;
            ady = dy > 0 ? dy : -dy;
            newEdge.C = (double)(s1.Coord.X * dx + s1.Coord.Y * dy + (dx * dx + dy * dy) * 0.5);

            if (adx > ady)
            {
                newEdge.A = 1.0;
                newEdge.B = dy / dx;
                newEdge.C /= dx;
            }
            else
            {
                newEdge.A = dx / dy;
                newEdge.B = 1.0;
                newEdge.C /= dy;
            }

            newEdge.EdgeNumber = nedges;
            nedges++;

            return newEdge;
        }

        // Assigns a unique site number to a given vertex and increments the total number of vertices.
        private void MakeVertex(Site v)
        {
            v.SiteNumber = nvertices;
            nvertices++;
        }

        // Initializes the priority queue by setting up the required data structures and initial values.
        private bool PQInitialize()
        {
            PQcount = 0;
            PQmin = 0;
            PQhashsize = 4 * sqrt_nsites;
            PQhash = new Halfedge[PQhashsize];

            for (int i = 0; i < PQhashsize; i++)
            {
                PQhash[i] = new Halfedge();
            }

            return true;
        }

        //  Calculates the bucket index in the priority queue for a given Halfedge based on its YStar value.
        private int PQBucket(Halfedge he)
        {
            int bucket;

            bucket = (int)((he.YStar - ymin) / deltay * PQhashsize);
            if (bucket < 0)
                bucket = 0;
            if (bucket >= PQhashsize)
                bucket = PQhashsize - 1;
            if (bucket < PQmin)
                PQmin = bucket;

            return bucket;
        }

        // Inserts a HalfEdge into the priority queue based on the given Site and offset.
        private void PQInsert(Halfedge he, Site v, double offset)
        {
            Halfedge last, next;

            he.Vertex = v;
            he.YStar = (double)(v.Coord.Y + offset);
            last = PQhash[PQBucket(he)];

            // Find the correct position in the linked list to insert the HalfEdge.
            while (
                (next = last.Next) != null
                &&
                (he.YStar > next.YStar || (he.YStar == next.YStar && v.Coord.X > next.Vertex.Coord.X))
            )
            {
                last = next;
            }

            // Insert the HalfEdge into the linked list.
            he.Next = last.Next;
            last.Next = he;
            PQcount++;
        }

        // Removes a HalfEdge from the priority queue.
        private void PQDelete(Halfedge he)
        {
            Halfedge last;

            if (he.Vertex != null)
            {
                last = PQhash[PQBucket(he)];
                // Find the HalfEdge to delete in the linked list.
                while (last.Next != he)
                {
                    last = last.Next;
                }

                // Remove the HalfEdge from the linked list.
                last.Next = he.Next;
                PQcount--;
                he.Vertex = null;
            }
        }

        private bool PQEmpty()
        {
            return (PQcount == 0);
        }

        private Point PQMin()
        {
            Point answer = new Point();

            while (PQhash[PQmin].Next == null)
            {
                PQmin++;
            }

            answer.X = PQhash[PQmin].Next.Vertex.Coord.X;
            answer.Y = PQhash[PQmin].Next.YStar;
            return answer;
        }

        private Halfedge PQExtractMin()
        {
            Halfedge curr;

            curr = PQhash[PQmin].Next;
            PQhash[PQmin].Next = curr.Next;
            PQcount--;

            return curr;
        }

        private Halfedge HECreate(Edge e, int pm)
        {
            Halfedge answer = new Halfedge();
            answer.Edge = e;
            answer.PM = pm;
            answer.Next = null;
            answer.Vertex = null;

            return answer;
        }


        private bool ELInitialize()
        {
            ELhashsize = 2 * sqrt_nsites;
            ELhash = new Halfedge[ELhashsize];

            for (int i = 0; i < ELhashsize; i++)
            {
                ELhash[i] = null;
            }

            ELleftend = HECreate(null, 0);
            ELrightend = HECreate(null, 0);
            ELleftend.Left = null;
            ELleftend.Right = ELrightend;
            ELrightend.Left = ELleftend;
            ELrightend.Right = null;
            ELhash[0] = ELleftend;
            ELhash[ELhashsize - 1] = ELrightend;

            return true;
        }

        private Halfedge ELRight(Halfedge he)
        {
            return he.Right;
        }

        private Halfedge ELLeft(Halfedge he)
        {
            return he.Left;
        }

        private Site LeftReg(Halfedge he)
        {
            if (he.Edge == null)
            {
                return bottomsite;
            }

            return (he.PM == LE ? he.Edge.Region1 : he.Edge.Region2);
        }

        private void ELInsert(Halfedge lb, Halfedge newHe)
        {
            newHe.Left = lb;
            newHe.Right = lb.Right;
            lb.Right.Left = newHe;
            lb.Right = newHe;
        }

        /* This delete routine can't reclaim the node, since pointers from the hash table
     * may be present.
     */
        private void ELDelete(Halfedge he)
        {
            he.Left.Right = he.Right;
            he.Right.Left = he.Left;
            he.Deleted = true;
        }

        // Get entry from hash table, pruning any deleted nodes
        private Halfedge ELGetHash(int b)
        {
            Halfedge he;
            if (b < 0 || b >= ELhashsize)
                return null;

            he = ELhash[b];
            if (he == null || !he.Deleted)
                return he;

            // Hash table points to deleted half edge. Patch as necessary.
            ELhash[b] = null;
            return null;
        }

        // Find the Halfedge that is immediately to the left of the input point
        private Halfedge FindHalfEdgeToLeftOfPoint(Point p)
        {
            int bucket;
            Halfedge he;

            // Use hash table to get close to desired halfedge
            bucket = (int)((p.X - xmin) / deltax * ELhashsize);

            // Make sure that the bucket position is within the range of the hash array
            if (bucket < 0) bucket = 0;
            if (bucket >= ELhashsize) bucket = ELhashsize - 1;

            he = ELGetHash(bucket);

            // If the HE isn't found, search backwards and forwards in the hash map
            // for the first non-null entry
            if (he == null)
            {
                for (int i = 1; i < ELhashsize; i++)
                {
                    if ((he = ELGetHash(bucket - i)) != null)
                        break;
                    if ((he = ELGetHash(bucket + i)) != null)
                        break;
                }
            }

            // Now search linear list of halfedges for the correct one
            if (he == ELleftend || (he != ELrightend && RightOf(he, p)))
            {
                // Keep going right on the list until either the end is reached, or
                // you find the 1st edge which the point isn't to the right of
                do
                {
                    he = he.Right;
                } while (he != ELrightend && RightOf(he, p));

                he = he.Left;
            }
            else
            {
                // If the point is to the left of the HalfEdge, then search left for
                // the HE just to the left of the point
                do
                {
                    he = he.Left;
                } while (he != ELleftend && !RightOf(he, p));
            }

            // Update hash table and reference counts
            if (bucket > 0 && bucket < ELhashsize - 1)
            {
                ELhash[bucket] = he;
            }

            return he;
        }


        private int k = 0;
        private void PushGraphEdge(Site leftSite, Site rightSite, double x1, double y1, double x2, double y2)
        {
            k++;
            if (k == 263) return;
            GraphEdge newEdge = new GraphEdge(x1, y1, x2, y2, leftSite.SiteNumber,  rightSite.SiteNumber);
            allEdges.Add(newEdge);
        }


        // Method used to clip edges of the Voronoi diagram to the borders of a given region
        private void ClipLine(Edge edge)
        {
            // Define variables to hold the borders of the region
            double regionLeftborder, regionRightBorder, regionLowerBorder, regionTopBorder;
            // Define variables to hold the sites that the edge was created from
            Site s1, s2;

            // Calculate the x and y coordinates of the two sites that the edge was created from
            double x1 = edge.Region1.Coord.X;
            double y1 = edge.Region1.Coord.Y;
            double x2 = edge.Region2.Coord.X;
            double y2 = edge.Region2.Coord.Y;
            double x = x2 - x1;
            double y = y2 - y1;

            // If the distance between the two sites is less than a certain value, ignore the edge
            if (Math.Sqrt((x * x) + (y * y)) < minDistanceBetweenSites)
            {
                return;
            }

            // Define the borders of the region
            regionLeftborder = borderMinX;
            regionLowerBorder = borderMinY;
            regionRightBorder = borderMaxX;
            regionTopBorder = borderMaxY;

            // Determine which site is on the left and which is on the right of the edge
            if (edge.A == 1.0 && edge.B >= 0.0)
            {
                s1 = edge.EndPoint2;
                s2 = edge.EndPoint1;
            }
            else
            {
                s1 = edge.EndPoint1;
                s2 = edge.EndPoint2;
            }

            // Handle the case where the edge is horizontal
            if (edge.A == 1.0)
            {
                // Set y1 to the lower border of the region
                y1 = regionLowerBorder;
                // If the first site is above the lower border, set y1 to its y coordinate
                if (s1 != null && s1.Coord.Y > regionLowerBorder)
                    y1 = s1.Coord.Y;
                // If y1 is greater than the upper border, set it to the upper border
                if (y1 > regionTopBorder)
                    y1 = regionTopBorder;
                // Calculate x1 based on the y1 value
                x1 = edge.C - edge.B * y1;
                // Set y2 to the upper border of the region
                y2 = regionTopBorder;
                // If the second site is below the upper border, set y2 to its y coordinate
                if (s2 != null && s2.Coord.Y < regionTopBorder)
                    y2 = s2.Coord.Y;
                // If y2 is less than the lower border, set it to the lower border
                if (y2 < regionLowerBorder)
                    y2 = regionLowerBorder;
                // Calculate x2 based on the y2 value
                x2 = edge.C - edge.B * y2;
                // If both x1 and x2 are outside the region, ignore the edge
                if (((x1 > regionRightBorder) & (x2 > regionRightBorder)) | ((x1 < regionLeftborder) & (x2 < regionLeftborder)))
                    return;

                // If x1 is greater than the right border, set it to the right border and recalculate y1
                if (x1 > regionRightBorder)
                {
                    x1 = regionRightBorder;
                    y1 = (edge.C - x1) / edge.B;
                }

                // If x1 is less than the left border, set it to the left border and recalculate y1
                if (x1 < regionLeftborder)
                {
                    x1 = regionLeftborder;
                    y1 = (edge.C - x1) / edge.B;
                }

                // If x2 is greater than the right border, set it to the right border and recalculate y2
                if (x2 > regionRightBorder)
                {
                    x2 = regionRightBorder;
                    y2 = (edge.C - x2) / edge.B;
                }

                // If x2 is less than the left border, set it to the left border and recalculate y2
                if (x2 < regionLeftborder)
                {
                    x2 = regionLeftborder;
                    y2 = (edge.C - x2) / edge.B;
                }
            }
            else // Handle the case where the edge is vertical
            {
                // Set x1 to the left border of the region
                x1 = regionLeftborder;
                // If the first site is to the right of the left border, set x1 to its x coordinate
                if (s1 != null && s1.Coord.X > regionLeftborder)
                    x1 = s1.Coord.X;
                // If x1 is greater than the right border, set it to the right border
                if (x1 > regionRightBorder)
                    x1 = regionRightBorder;
                // Calculate y1 based on the x1 value
                y1 = edge.C - edge.A * x1;

                // Set x2 to the right border of the region
                x2 = regionRightBorder;
                // If the second site is to the left of the right border, set x2 to its x coordinate
                if (s2 != null && s2.Coord.X < regionRightBorder)
                    x2 = s2.Coord.X;
                // If x2 is less than the left border, set it to the left border
                if (x2 < regionLeftborder)
                    x2 = regionLeftborder;
                // Calculate y2 based on the x2 value
                y2 = edge.C - edge.A * x2;

                // If both y1 and y2 are outside the region, ignore the edge
                if (((y1 > regionTopBorder) & (y2 > regionTopBorder)) | ((y1 < regionLowerBorder) & (y2 < regionLowerBorder)))
                    return;

                // If y1 is greater than the upper border, set it to the upper border and recalculate x1
                if (y1 > regionTopBorder)
                {
                    y1 = regionTopBorder;
                    x1 = (edge.C - y1) / edge.A;
                }

                // If y1 is less than the lower border, set it to the lower border and recalculate x1
                if (y1 < regionLowerBorder)
                {
                    y1 = regionLowerBorder;
                    x1 = (edge.C - y1) / edge.A;
                }

                // If y2 is greater than the upper border, set it to the upper border and recalculate x2
                if (y2 > regionTopBorder)
                {
                    y2 = regionTopBorder;
                    x2 = (edge.C - y2) / edge.A;
                }

                // If y2 is less than the lower border, set it to the lower border and recalculate x2
                if (y2 < regionLowerBorder)
                {
                    y2 = regionLowerBorder;
                    x2 = (edge.C - y2) / edge.A;
                }
            }

            // Add the clipped edge to the Voronoi diagram
            PushGraphEdge(edge.Region1, edge.Region2, x1, y1, x2, y2);
        }

        private void SetEndpoint(Edge edge, int leftRight, Site site)
        {
            if (leftRight == LE)
            {
                if (edge.EndPoint1 == null)
                {
                    edge.EndPoint1 = site;
                }
                else
                {
                    edge.EndPoint2 = site;
                }
            }
            else
            {
                if (edge.EndPoint2 == null)
                {
                    edge.EndPoint2 = site;
                }
                else
                {
                    edge.EndPoint1 = site;
                }
            }

            if (edge.EndPoint1 != null && edge.EndPoint2 != null)
            {
                ClipLine(edge);
            }
        }

        private Site RightSite(Halfedge halfedge)
        {
            if (halfedge.Edge == null)
            {
                // If this halfedge has no edge, return the bottom site
                return bottomsite;
            }

            // If the PM field is LE, return the site 1 that this edge bisects,
            // otherwise return site 0
            return (halfedge.PM == LE ? halfedge.Edge.Region2 : halfedge.Edge.Region1);
        }

        private bool RightOf(Halfedge halfedge, Point point)
        {
            Edge edge = halfedge.Edge;
            Site topSite = edge.Region2;
            bool rightOfSite;

            if (point.X > topSite.Coord.X)
                rightOfSite = true;
            else
                rightOfSite = false;

            if (rightOfSite && halfedge.PM == LE)
                return true;
            if (!rightOfSite && halfedge.PM == RE)
                return false;

            bool isAbove, fastCheck;
            double pointDx, pointDy, siteDx, temp1, temp2, temp3, yl;

            if (edge.A == 1.0)
            {
                pointDx = point.X - topSite.Coord.X;
                pointDy = point.Y - topSite.Coord.Y;
                fastCheck = false;

                if ((!rightOfSite & (edge.B < 0.0)) || (rightOfSite & (edge.B >= 0.0)))
                {
                    isAbove = pointDy >= edge.B * pointDx;
                    fastCheck = isAbove;
                }
                else
                {
                    isAbove = point.X + point.Y * edge.B > edge.C;
                    if (edge.B < 0.0)
                        isAbove = !isAbove;
                    if (!isAbove)
                        fastCheck = true;
                }

                if (!fastCheck)
                {
                    siteDx = topSite.Coord.X - edge.Region1.Coord.X;
                    isAbove = edge.B * (pointDx * pointDx - pointDy * pointDy)
                              < siteDx * pointDy * (1.0 + 2.0 * pointDx / siteDx + edge.B * edge.B);

                    if (edge.B < 0)
                        isAbove = !isAbove;
                }
            }
            else // edge.B == 1.0
            {
                yl = edge.C - edge.A * point.X;
                temp1 = point.Y - yl;
                temp2 = point.X - topSite.Coord.X;
                temp3 = yl - topSite.Coord.Y;
                isAbove = temp1 * temp1 > temp2 * temp2 + temp3 * temp3;
            }

            return (halfedge.PM == LE ? isAbove : !isAbove);
        }

        private Site rightreg(Halfedge he)
        {
            if (he.Edge == null)
                // if this halfedge has no edge, return the bottom site (whatever
                // that is)
            {
                return (bottomsite);
            }

            // if the ELpm field is zero, return the site 0 that this edge bisects,
            // otherwise return site number 1
            return (he.PM == LE ? he.Edge.Region2 : he.Edge.Region1);
        }

        private double dist(Site s, Site t)
        {
            double dx, dy;
            dx = s.Coord.X - t.Coord.X;
            dy = s.Coord.Y - t.Coord.Y;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        // create a new site where the HalfEdges el1 and el2 intersect - note that
        // the Point in the argument list is not used, don't know why it's there
        private Site intersect(Halfedge el1, Halfedge el2)
        {
            Edge e1, e2, e;
            Halfedge el;
            double d, xint, yint;
            bool right_of_site;
            Site v; // vertex

            e1 = el1.Edge;
            e2 = el2.Edge;

            if (e1 == null || e2 == null)
                return null;

            // if the two edges bisect the same parent, return null
            if (e1.Region2 == e2.Region2)
                return null;

            d = e1.A * e2.B - e1.B * e2.A;
            if (-1.0e-10 < d && d < 1.0e-10)
                return null;

            xint = (e1.C * e2.B - e2.C * e1.B) / d;
            yint = (e2.C * e1.A - e1.C * e2.A) / d;

            if ((e1.Region2.Coord.Y < e2.Region2.Coord.Y)
                || (e1.Region2.Coord.Y == e2.Region2.Coord.Y && e1.Region2.Coord.X < e2.Region2.Coord.X))
            {
                el = el1;
                e = e1;
            }
            else
            {
                el = el2;
                e = e2;
            }

            right_of_site = xint >= e.Region2.Coord.X;
            if ((right_of_site && el.PM == LE)
                || (!right_of_site && el.PM == RE))
                return null;

            // create a new site at the point of intersection - this is a new vector
            // event waiting to happen
            v = new Site();
            v.Coord.X = xint;
            v.Coord.Y = yint;
            return v;
        }

        /*
         * implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax, deltax,
         * deltay (can all be estimates). Performance suffers if they are wrong;
         * better to make nsites, deltax, and deltay too big than too small. (?)
         */
        private bool VoronoiBd()
        {
            Site newsite, bot, top, temp, p;
            Site v;
            Point newintstar = null;
            int pm;
            Halfedge lbnd, rbnd, llbnd, rrbnd, bisector;
            Edge e;

            PQInitialize();
            ELInitialize();

            bottomsite = NextOne();
            newsite = NextOne();
            while (true)
            {
                if (!PQEmpty())
                {
                    newintstar = PQMin();
                }
                // if the lowest site has a smaller y value than the lowest vector
                // intersection,
                // process the site otherwise process the vector intersection

                if (newsite != null && (PQEmpty()
                                        || newsite.Coord.Y < newintstar.Y
                                        || (newsite.Coord.Y == newintstar.Y
                                            && newsite.Coord.X < newintstar.X)))
                {
                    /* new site is smallest -this is a site event */
                    // get the first HalfEdge to the LEFT of the new site
                    lbnd = FindHalfEdgeToLeftOfPoint((newsite.Coord));
                    // get the first HalfEdge to the RIGHT of the new site
                    rbnd = ELRight(lbnd);
                    // if this halfedge has no edge,bot =bottom site (whatever that
                    // is)
                    bot = rightreg(lbnd);
                    // create a new edge that bisects
                    e = Bisect(bot, newsite);

                    // create a new HalfEdge, setting its ELpm field to 0
                    bisector = HECreate(e, LE);
                    // insert this new bisector edge between the left and right
                    // vectors in a linked list
                    ELInsert(lbnd, bisector);

                    // if the new bisector intersects with the left edge,
                    // remove the left edge's vertex, and put in the new one
                    if ((p = intersect(lbnd, bisector)) != null)
                    {
                        PQDelete(lbnd);
                        PQInsert(lbnd, p, dist(p, newsite));
                    }

                    lbnd = bisector;
                    // create a new HalfEdge, setting its ELpm field to 1
                    bisector = HECreate(e, RE);
                    // insert the new HE to the right of the original bisector
                    // earlier in the IF stmt
                    ELInsert(lbnd, bisector);

                    // if this new bisector intersects with the new HalfEdge
                    if ((p = intersect(bisector, rbnd)) != null)
                    {
                        // push the HE into the ordered linked list of vertices
                        PQInsert(bisector, p, dist(p, newsite));
                    }

                    newsite = NextOne();
                }
                else if (!PQEmpty())
                    /* intersection is smallest - this is a vector event */
                {
                    // pop the HalfEdge with the lowest vector off the ordered list
                    // of vectors
                    lbnd = PQExtractMin();
                    // get the HalfEdge to the left of the above HE
                    llbnd = ELLeft(lbnd);
                    // get the HalfEdge to the right of the above HE
                    rbnd = ELRight(lbnd);
                    // get the HalfEdge to the right of the HE to the right of the
                    // lowest HE
                    rrbnd = ELRight(rbnd);
                    // get the Site to the left of the left HE which it bisects
                    bot = LeftReg(lbnd);
                    // get the Site to the right of the right HE which it bisects
                    top = rightreg(rbnd);

                    v = lbnd.Vertex; // get the vertex that caused this event
                    MakeVertex(v); // set the vertex number - couldn't do this
                    // earlier since we didn't know when it would be processed
                    SetEndpoint(lbnd.Edge, lbnd.PM, v);
                    // set the endpoint of
                    // the left HalfEdge to be this vector
                    SetEndpoint(rbnd.Edge, rbnd.PM, v);
                    // set the endpoint of the right HalfEdge to
                    // be this vector
                    ELDelete(lbnd); // mark the lowest HE for
                    // deletion - can't delete yet because there might be pointers
                    // to it in Hash Map
                    PQDelete(rbnd);
                    // remove all vertex events to do with the right HE
                    ELDelete(rbnd); // mark the right HE for
                    // deletion - can't delete yet because there might be pointers
                    // to it in Hash Map
                    pm = LE; // set the pm variable to zero

                    if (bot.Coord.Y > top.Coord.Y)
                        // if the site to the left of the event is higher than the
                        // Site
                    {
                        // to the right of it, then swap them and set the 'pm'
                        // variable to 1
                        temp = bot;
                        bot = top;
                        top = temp;
                        pm = RE;
                    }

                    e = Bisect(bot, top); // create an Edge (or line)
                    // that is between the two Sites. This creates the formula of
                    // the line, and assigns a line number to it
                    bisector = HECreate(e, pm); // create a HE from the Edge 'e',
                    // and make it point to that edge
                    // with its ELedge field
                    ELInsert(llbnd, bisector); // insert the new bisector to the
                    // right of the left HE
                    SetEndpoint(e, RE - pm, v); // set one endpoint to the new edge
                    // to be the vector point 'v'.
                    // If the site to the left of this bisector is higher than the
                    // right Site, then this endpoint
                    // is put in position 0; otherwise in pos 1

                    // if left HE and the new bisector intersect, then delete
                    // the left HE, and reinsert it
                    if ((p = intersect(llbnd, bisector)) != null)
                    {
                        PQDelete(llbnd);
                        PQInsert(llbnd, p, dist(p, bot));
                    }

                    // if right HE and the new bisector intersect, then
                    // reinsert it
                    if ((p = intersect(bisector, rrbnd)) != null)
                    {
                        PQInsert(bisector, p, dist(p, bot));
                    }
                }
                else
                {
                    break;
                }
            }

            for (lbnd = ELRight(ELleftend); lbnd != ELrightend; lbnd = ELRight(lbnd))
            {
                e = lbnd.Edge;
                ClipLine(e);
            }

            return true;
        }
    } // Voronoi Class End
} // namespace Voronoi2 End