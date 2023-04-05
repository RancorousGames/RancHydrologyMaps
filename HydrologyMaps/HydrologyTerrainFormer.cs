using System.Numerics;

namespace HydrologyMaps;

public class HydrologyTerrainFormer
{
    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float riverDepth,
        float baseRiverWidth, float falloffFactor, Random random)
    {
        foreach (var root in riverRoots)
        {
            CarveRiverDFS(heightMap, root, riverDepth, baseRiverWidth, falloffFactor, random);
        }
    }

    private static void CarveRiverDFS(float[,] heightMap, DirectedNode currentNode, float riverDepth,
        float baseRiverWidth, float falloffFactor, Random random)
    {
        if (currentNode.Children.Count == 0) return;

        foreach (var childGraphNode in currentNode.Children)
        {
            var childNode = (DirectedNode)childGraphNode;
            double parentFlowRate = currentNode.FlowRate / 25;
            double childFlowRate = childNode.FlowRate / 25;
            float parentWidth = baseRiverWidth * (float)Math.Pow(parentFlowRate, 0.3);
            float childWidth = baseRiverWidth * (float)Math.Pow(childFlowRate, 0.3);

            CarveRiverSegment(heightMap, currentNode.Point, childNode.Point, riverDepth, parentWidth, childWidth,
                falloffFactor, random);
            CarveRiverDFS(heightMap, (DirectedNode)childNode, riverDepth, baseRiverWidth, falloffFactor, random);
        }
    }

    private static void CarveRiverSegment(float[,] heightMap, Point start, Point end, float depth, float startWidth,
        float endWidth, float falloffFactor, Random random)
    {
        // Calculate the line equation between start and end points
        Vector2 direction = (end - start).ToVector();
        int numSteps = (int)direction.Length();
        direction /= numSteps;

        // Iterate over the points along the line between start and end points
        for (int i = 0; i <= numSteps; i++)
        {
            Point position = start + Point.FromVector(direction * i);
            int x = Math.Clamp(position.X, 0, heightMap.GetLength(0) - 1);
            int y = Math.Clamp(position.Y, 0, heightMap.GetLength(1) - 1);

            float t = (float)i / numSteps;
            float baseWidth = Lerp(startWidth, endWidth, t);

            // Apply Perlin noise to the river width
            float perlinSeed = random.Next(10000);
            float perlinValue = 0.1f; // (float)PerlinNoise.GetValue(position.X / 100f + perlinSeed, position.Y / 100f + perlinSeed, 0);
            float currentWidth = baseWidth * (0.8f + 0.2f * perlinValue);

            // Update the heightMap value by subtracting the depth and considering the width
            for (int offsetX = -((int)currentWidth); offsetX <= (int)currentWidth; offsetX++)
            {
                for (int offsetY = -((int)currentWidth); offsetY <= (int)currentWidth; offsetY++)
                {
                    int nx = Math.Clamp(x + offsetX, 0, heightMap.GetLength(0) - 1);
                    int ny = Math.Clamp(y + offsetY, 0, heightMap.GetLength(1) - 1);
                    float distance = (float)Math.Sqrt(offsetX * offsetX + offsetY * offsetY);


                    if (distance <= currentWidth)
                    {
                        // Calculate depth falloff using the falloffFactor
                        float depthFalloff = 1- (float)Math.Pow(distance / currentWidth, falloffFactor);
                        float currentDepth = depth * depthFalloff;

                        heightMap[nx, ny] = Math.Max(heightMap[nx, ny] - currentDepth, 0);
                    }
                }
            }
        }
    }


    /*
    // Carve rivers with smoothed paths, width variation, and curvature
    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float riverDepth,
        float falloff, float smoothness, float minWidth, float maxWidth)
    {
        if (smoothness <= 0) smoothness = 0.1f;

        foreach (var root in riverRoots)
        {
            CarveRiverDFS(heightMap, root, riverDepth, falloff, smoothness, minWidth, maxWidth);
        }
    }

    private static void CarveRiverDFS(float[,] heightMap, DirectedNode currentNode, float riverDepth,
        float falloff, float smoothness, float minWidth, float maxWidth)
    {
        if (currentNode.Children.Count == 0) return;

        foreach (var childNode in currentNode.Children)
        {
            List<DirectedNode> fullPath = BuildFullPath(currentNode);
            List<Point> smoothedPath = SmoothRiverPath(fullPath, smoothness);

            double avgFlowRate = Lerp(currentNode.FlowRate, ((DirectedNode)childNode).FlowRate, 0.5f);

            for (int i = 0; i < smoothedPath.Count - 1; i++)
            {
                Point start = smoothedPath[i];
                Point end = smoothedPath[i + 1];

                double depth = riverDepth * avgFlowRate;
                double width = Lerp(minWidth, maxWidth, avgFlowRate);

                CarveRiverSegment(heightMap, start, end, (float)depth, (float)width, falloff);
            }

            // Recursively carve the river segments of child nodes
            CarveRiverDFS(heightMap, (DirectedNode)childNode, riverDepth, falloff, smoothness, minWidth, maxWidth);
        }
    }

    private static List<DirectedNode> BuildFullPath(DirectedNode currentNode)
    {
        List<DirectedNode> fullPath = new List<DirectedNode> { currentNode };

        DirectedNode parentNode = (DirectedNode)currentNode.Parent;
        while (parentNode != null)
        {
            fullPath.Insert(0, parentNode);
            parentNode = (DirectedNode)parentNode.Parent;
        }

        return fullPath;
    }

    // Smooth the river path using Catmull-Rom splines
    private static List<Point> SmoothRiverPath(List<DirectedNode> path, float smoothness)
    {
        List<Point> smoothedPath = new List<Point>();
        for (float t = 0; t < path.Count - 1; t += smoothness)
        {
            int index = (int)Math.Floor(t);
            Vector2 point = CatmullRomSpline(path[index].PositionVector, path[index + 1].PositionVector,
                path[Math.Clamp(index + 2, 0, path.Count - 1)].PositionVector,
                path[Math.Clamp(index - 1, 0, path.Count - 1)].PositionVector, t - index);
            smoothedPath.Add(new Point((int)point.X, (int)point.Y));
        }

        smoothedPath.Add(path[^1].Point);
        return smoothedPath;
    }

    // Calculate the Catmull-Rom spline for a set of control points and a parameter t
    private static Vector2 CatmullRomSpline(Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3, float t)
    {
        Vector2 a = 2f * p1;
        Vector2 b = p2 - p0;
        Vector2 c = 2f * p0 - 5f * p1 + 4f * p2 - p3;
        Vector2 d = -p0 + 3f * p1 - 3f *
            p2 + p3;

        return 0.5f * (a + (b * t) + (c * t * t) + (d * t * t * t));
    }

    // Carve a river segment between start and end points, using depth, width, and falloff
    private static void CarveRiverSegment(float[,] heightMap, Point start, Point end, float depth,
        float width, float falloff)
    {
        Point dirPoint = (end - start);
        Vector2 direction = new Vector2(dirPoint.X, dirPoint.Y);
        direction /= direction.Length();

        Vector2 side = new Vector2(-direction.Y, direction.X);

        int widthInt = (int)Math.Ceiling(width);

        for (float t = 0; t <= 1; t += 1 / falloff)
        {
            Point center = Point.Lerp(start, end, t);
            for (int y = -widthInt; y <= widthInt; y++)
            {
                for (int x = -widthInt; x <= widthInt; x++)
                {
                    Vector2 offset = x * direction + y * side;
                    int posX = center.X + (int)(offset.X);
                    int posY = center.Y + (int)(offset.Y);

                    posX = Math.Clamp(posX, 0, heightMap.GetLength(0) - 1);
                    posY = Math.Clamp(posY, 0, heightMap.GetLength(1) - 1);

                    float distance = new Vector2(x, y).Length();
                    if (distance <= falloff)
                    {
                        float riverElevation = (float)Lerp(depth, 0, distance / falloff);
                        heightMap[posX, posY] = Math.Min(heightMap[posX, posY], riverElevation);
                    }
                }
            }
        }
    }
*/
    public static float Lerp(float a, float b, float t)
    {
        return a + (b - a) * t;
    }
}