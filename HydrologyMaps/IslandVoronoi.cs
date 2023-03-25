namespace HydrologyMaps;

using System.Numerics;
using Voronoi2;

public class IslandVoronoi
{
    public static List<GraphEdge> GenerateVoronoiEdges(List<DirectedNode> riverNodes)
    {
        var voroObject = new Voronoi(0.1);

        //double[] xVal = riverNodes.Select(x => (double)x.X).ToArray();
        //double[] yVal = riverNodes.Select(x => (double)x.Y).ToArray();

        List<GraphEdge> voronoi = voroObject.GenerateVoronoi(riverNodes, 0, 512, 0, 512);

        return voronoi;
    }

    public static List<GraphEdge> GenerateExtraEdges(List<DirectedNode> borderNodes, List<GraphEdge> voronoiEdges,
        float[,] heightMap)
    {
        // Step 1: Connect border nodes to form boundary edges
        List<GraphEdge> boundaryEdges = new List<GraphEdge>();
        for (int i = 0; i < borderNodes.Count; i++)
        {
            Point start = borderNodes[i].Point;
            Point end = borderNodes[(i + 1) % borderNodes.Count].Point;
            boundaryEdges.Add(new GraphEdge(start.X, start.Y, end.X, end.Y));
        }

        int k = 0;

        var firstRiverNodeX = boundaryEdges[0].X1;
        var firstRiverNodeY = boundaryEdges[0].Y1;
        
        
        // remove all voronoi edges that are entirely in the ocean
        voronoiEdges = voronoiEdges
            .Where(e => heightMap[(int)e.X1, (int)e.Y1] > 0 || heightMap[(int)e.X2, (int)e.Y2] > 0).ToList();
        List<GraphEdge> allEdges = new List<GraphEdge>(voronoiEdges);
        
        foreach (GraphEdge boundaryEdge in boundaryEdges)
        {
            List<GraphEdge> newEdges = new List<GraphEdge>();

            bool intersectionFound = false;
            foreach (GraphEdge voronoiEdge in voronoiEdges)
            {
                
                
                Voronoi2.Point? intersection = boundaryEdge.Intersection(voronoiEdge);
                if (intersection != null)
                {
                    // We have a voronoi edge that crosses our clockwise direction going from oceanmouth to oceanmouth, now cut it in half and pick the inland one to keep
                    Vector2 direction = new Vector2((float)(boundaryEdge.X2 - boundaryEdge.X1),
                        (float)(boundaryEdge.Y2 - boundaryEdge.Y1));
                    var toStartSiteX = voronoiEdge.X1 - boundaryEdge.X1;
                    var toStartSiteY = voronoiEdge.Y1 - boundaryEdge.Y1;
                    var crossProduct = (direction.X * toStartSiteY) - (toStartSiteX * direction.Y);

                    if (crossProduct > 3)
                        // If the start of the island-border-crossing voronoi edge is to the left of our clockwise direction then add the half of it that uses... start?
                        allEdges.Add(new GraphEdge(voronoiEdge.X1, voronoiEdge.Y1, intersection.X, intersection.Y));
                    else
                        // else add the other half 
                        allEdges.Add(new GraphEdge(voronoiEdge.X2, voronoiEdge.Y2, intersection.X, intersection.Y));

                    // lets keep the oceanmouth to oceanmouth edge (split in two which might be important later, not sure)
                    allEdges.Add(new GraphEdge(boundaryEdge.X1, boundaryEdge.Y1, intersection.X, intersection.Y));
                    allEdges.Add(new GraphEdge(boundaryEdge.X2, boundaryEdge.Y2, intersection.X, intersection.Y));

                    // remove the full long edge now that we added the half
                    allEdges.Remove(voronoiEdge);
                    if (boundaryEdge.X2 == firstRiverNodeX && boundaryEdge.Y2 == firstRiverNodeY)
                    {
                        return allEdges;
                    }

                    intersectionFound = true;
                    break;
                }

            }
                
            // If no intersection, lets connect the boundary edge anyway
            if (!intersectionFound) allEdges.Add(boundaryEdge);
        }

        // Step 5: Return set of filtered edges
        return allEdges;
    }

    private static bool IsSiteConnectedToSite((double x, double y) site1, (double x, double y) site2,
        List<GraphEdge> edges)
    {
        // Check if there is an edge that connects the two sites
        foreach (GraphEdge edge in edges)
        {
            if ((edge.X1.Equals(site1.x) && edge.Y1.Equals(site1.y) && edge.X2.Equals(site2.x) &&
                 edge.Y2.Equals(site2.y)) ||
                (edge.X1.Equals(site2.x) && edge.Y1.Equals(site2.y) && edge.X2.Equals(site1.x) &&
                 edge.Y2.Equals(site1.y)))
            {
                return true;
            }
        }

        // If no connecting edge is found, return false
        return false;
    }
}