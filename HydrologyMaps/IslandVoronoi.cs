namespace HydrologyMaps;
using System.Numerics;
using Voronoi2;
using Point = System.Drawing.Point;

public class IslandVoronoi
{
    public static List<Edge> GenerateExtraEdges(List<DirectedNode> borderNodes, List<Edge> voronoiEdges)
{
    // Step 1: Connect border nodes to form boundary edges
    List<Edge> boundaryEdges = new List<Edge>();
    for (int i = 0; i < borderNodes.Count; i++)
    {
        Point start = borderNodes[i].Point;
        Point end = borderNodes[(i + 1) % borderNodes.Count].Point;
        boundaryEdges.Add(new Edge(new Site(start.X, start.Y), new Site(end.X, end.Y)));
    }

    // Step 2: Connect last and first border nodes if not already connected
    Site firstSite = boundaryEdges[0].LeftSite;
    Site lastSite = boundaryEdges[boundaryEdges.Count - 1].RightSite;
    if (!IsSiteConnectedToSite(firstSite, lastSite, voronoiEdges))
    {
        boundaryEdges.Add(new Edge(firstSite, lastSite));
    }

    // Step 3: Split Voronoi edges at boundary intersections
    List<Edge> allEdges = new List<Edge>();
    foreach (Edge voronoiEdge in voronoiEdges)
    {
        List<Edge> newEdges = new List<Edge>();
        foreach (Edge boundaryEdge in boundaryEdges)
        {
            Point intersection = boundaryEdge.Intersection(voronoiEdge);
            if (intersection != null)
            {
                Site startSite = voronoiEdge.LeftSite;
                Site endSite = voronoiEdge.RightSite;
                if (voronoiEdge.IsPointToLeft(intersection))
                {
                    startSite = voronoiEdge.RightSite;
                    endSite = voronoiEdge.LeftSite;
                }

                Vector2 direction = new Vector2(intersection.X - startSite.X, intersection.Y - startSite.Y);
                direction.Normalize();

                Edge newEdge1 = new Edge(startSite, new Site(intersection.X, intersection.Y));
                newEdge1.Direction = direction;
                Edge newEdge2 = new Edge(endSite, new Site(intersection.X, intersection.Y));
                newEdge2.Direction = -direction;

                newEdges.Add(newEdge1);
                newEdges.Add(newEdge2);
            }
        }

        if (newEdges.Count > 0)
        {
            // Add new edges if Voronoi edge intersects a boundary edge
            allEdges.AddRange(newEdges);
        }
        else
        {
            // Add original Voronoi edge if no intersection
            allEdges.Add(voronoiEdge);
        }
    }

    // Step 4: Return set of all edges
    return allEdges;
}

private static bool IsSiteConnectedToSite(Site site1, Site site2, List<Edge> edges)
{
    foreach (Edge edge in edges)
    {
        if ((edge.LeftSite == site1 && edge.RightSite == site2) ||
}