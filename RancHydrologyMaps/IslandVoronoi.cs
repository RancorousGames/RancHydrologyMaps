namespace HydrologyMaps;

using Voronoi2;

public class IslandVoronoi
{
    public static List<GraphEdge> GenerateVoronoiEdges(List<IGraphNode> riverNodes)
    {
        var voroObject = new Voronoi(0.1);

        List<GraphEdge> voronoi = voroObject.GenerateVoronoi(riverNodes, 0, 511, 0, 511);

        return voronoi;
    }

    public static List<GraphEdge> GenerateExtraEdges(List<DirectedNode> borderNodes, List<GraphEdge> voronoiEdges,
        float[,] heightMap)
    {
        // Step 1: Connect border nodes to form boundary edges
        List<GraphEdge> contourEdges = new List<GraphEdge>();
        for (int i = 0; i < borderNodes.Count; i++)
        {
            Point start = borderNodes[i].Point;
            Point end = borderNodes[(i + 1) % borderNodes.Count].Point;
            contourEdges.Add(new GraphEdge(start.X, start.Y, end.X, end.Y));
        }
        
        var firstRiverNodeX = contourEdges[0].X1;
        var firstRiverNodeY = contourEdges[0].Y1;

        // remove all voronoi edges that are entirely in the ocean
        voronoiEdges = voronoiEdges
            .Where(e => heightMap[(int)e.X1, (int)e.Y1] > 0.15 || heightMap[(int)e.X2, (int)e.Y2] > 0.15).ToList();

        List<GraphEdge> internalEdges = new List<GraphEdge>(voronoiEdges);
        List<GraphEdge> contourEdgesToAdd = new List<GraphEdge>();
        foreach (GraphEdge boundaryEdge in contourEdges)
        {
            List<GraphEdge> newEdges = new List<GraphEdge>();

            bool intersectionFound = false;
            foreach (GraphEdge voronoiEdge in voronoiEdges)
            {
                Voronoi2.Point? intersection = boundaryEdge.Intersection(voronoiEdge);
                if (intersection != null)
                {
                    // We have a voronoi edge that crosses our clockwise direction going from oceanmouth to oceanmouth, now cut it in half and pick the inland one to keep
                    Vector2D direction = new Vector2D(
                        (boundaryEdge.X2 - boundaryEdge.X1),
                        (boundaryEdge.Y2 - boundaryEdge.Y1));
                    var toStartSiteX = voronoiEdge.X1 - boundaryEdge.X1;
                    var toStartSiteY = voronoiEdge.Y1 - boundaryEdge.Y1;
                    var crossProduct = (direction.X * toStartSiteY) - (toStartSiteX * direction.Y);

                    var e = voronoiEdge;
                    if (crossProduct > 0)
                        // If the start of the island-border-crossing voronoi edge is to the left of our clockwise direction then add the half of it that uses... start?
                        internalEdges.Add(new GraphEdge(voronoiEdge.X1, voronoiEdge.Y1, intersection.X,
                            intersection.Y));
                    else if (crossProduct < 0)
                        // else add the other half 
                        internalEdges.Add(new GraphEdge(voronoiEdge.X2, voronoiEdge.Y2, intersection.X,
                            intersection.Y));

                    // lets keep the contouredges edge (split in two which might be important later, not sure)
                    contourEdgesToAdd.Add(new GraphEdge(boundaryEdge.X1, boundaryEdge.Y1, intersection.X,
                        intersection.Y));
                    contourEdgesToAdd.Add(new GraphEdge(boundaryEdge.X2, boundaryEdge.Y2, intersection.X,
                        intersection.Y));

                    // remove the full long edge now that we added the half
                    internalEdges.Remove(voronoiEdge);

                    intersectionFound = true;
                }
            }
            
            // ReSharper disable CompareOfFloatsByEqualityOperator
            if (boundaryEdge.X2 == firstRiverNodeX && boundaryEdge.Y2 == firstRiverNodeY)
            {
                break;
            }
            // ReSharper restore CompareOfFloatsByEqualityOperator

            // If no intersection, lets connect the boundary edge anyway
            if (!intersectionFound) contourEdgesToAdd.Add(boundaryEdge);
        }

        // Step 5: Return set of filtered edges

        return FilterValidEdgesAndAdd(heightMap, internalEdges, contourEdges, contourEdgesToAdd);
    }

    private static List<GraphEdge> FilterValidEdgesAndAdd(float[,] heightMap, List<GraphEdge> edgesToFilter,
        List<GraphEdge> contourEdges, List<GraphEdge> edgesToAdd)
    {
        List<GraphEdge> filteredEdges = new List<GraphEdge>(edgesToFilter.Count);

        foreach (GraphEdge e in edgesToFilter)
        {
            if (heightMap[(int)e.X1, (int)e.Y1] > 0 && heightMap[(int)e.X2, (int)e.Y2] > 0)
            {
                filteredEdges.Add(e);
            }
            else if (contourEdges.Any(ce => ce.Intersection(e) != null))
            {
                filteredEdges.Add(e);
            }
        }

        edgesToAdd.ForEach(e => filteredEdges.Add(e));

        return filteredEdges;
    }
}