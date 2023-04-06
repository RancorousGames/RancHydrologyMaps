using System.Numerics;

namespace HydrologyMaps;

public class HydrologyTerrainFormer
{
    private static FastNoiseLite _noiseGenerator = new();

    public HydrologyTerrainFormer()
    {
        _noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        _noiseGenerator.SetFrequency(0.003f);
        _noiseGenerator.SetFractalGain(0);
    }

    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float baseRiverWidth,
        float riverDepth, float smoothness, float falloff)
    {
        float[,] accumulatedDepthMap = new float[heightMap.GetLength(0), heightMap.GetLength(1)];

        foreach (var root in riverRoots)
        {
            AccumulateRiverSegmentsDFS(accumulatedDepthMap, root, baseRiverWidth, riverDepth, smoothness, falloff);
        }

        ApplyAccumulatedDepthMap(heightMap, accumulatedDepthMap);
    }

    private static void AccumulateRiverSegmentsDFS(float[,] accumulatedDepthMap, DirectedNode currentNode,
        float baseRiverWidth, float riverDepth, float smoothness, float falloff)
    {
        foreach (var childGraphNode in currentNode.Children)
        {
            var childNode = (DirectedNode)childGraphNode;

            List<Point> fullPath = BuildFullPath(currentNode);
            List<Point> smoothedPath = GenerateSmoothedPath(fullPath, smoothness);

            float startFlowRate = (float)currentNode.FlowRate / 25;
            float endFlowRate = (float)childNode.FlowRate / 25;

            AccumulateRiverSegment(accumulatedDepthMap, smoothedPath, startFlowRate, endFlowRate, baseRiverWidth,
                riverDepth, falloff);

            AccumulateRiverSegmentsDFS(accumulatedDepthMap, childNode, baseRiverWidth, riverDepth, smoothness, falloff);
        }

        // If the current node is a leaf node, process it here.
        if (currentNode.Children.Count == 0)
        {
            ProcessLeafNode(accumulatedDepthMap, currentNode, baseRiverWidth, riverDepth, smoothness, falloff);
        }
    }

    private static void ProcessLeafNode(float[,] accumulatedDepthMap, DirectedNode leafNode, float baseRiverWidth,
        float riverDepth, float smoothness, float falloff)
    {
        if (leafNode.Parent == null) return;

        DirectedNode parentNode = (DirectedNode)leafNode.Parent;

        List<Point> fullPath = new List<Point> { parentNode.Point, leafNode.Point };
        List<Point> smoothedPath = GenerateSmoothedPath(fullPath, smoothness);

        float startFlowRate = (float)parentNode.FlowRate / 25;
        float endFlowRate = (float)leafNode.FlowRate / 25;

        AccumulateRiverSegment(accumulatedDepthMap, smoothedPath, startFlowRate, endFlowRate, baseRiverWidth,
            riverDepth, falloff);
    }

    private static List<Point> BuildFullPath(DirectedNode currentNode)
    {
        List<DirectedNode> fullPath = new List<DirectedNode> { currentNode };

        DirectedNode parentNode = (DirectedNode)currentNode.Parent;
        while (parentNode != null)
        {
            fullPath.Insert(0, parentNode);
            parentNode = (DirectedNode)parentNode.Parent;
        }

        return fullPath.Select(node => node.Point).ToList();
    }

    private static List<Point> GenerateSmoothedPath(List<Point> fullPath, float smoothness)
    {
        // todo: Catmull-Rom splines
        return fullPath;
    }

    private static void AccumulateRiverSegment(float[,] accumulatedDepthMap, List<Point> smoothedPath,
        float startFlowRate, float endFlowRate, float baseRiverWidth, float riverDepth, float falloffFactor)
    {
        for (int i = 0; i < smoothedPath.Count - 1; i++)
        {
            Point start = smoothedPath[i];
            Point end = smoothedPath[i + 1];

            float t = (float)i / (smoothedPath.Count - 1);

            float startWidth = baseRiverWidth * (float)Math.Pow(startFlowRate, 0.3);
            float endWidth = baseRiverWidth * (float)Math.Pow(endFlowRate, 0.3);

            float startDepth = riverDepth * startFlowRate;
            float endDepth = riverDepth * endFlowRate;
            
            CarveRiverSegment(accumulatedDepthMap, start, end, startDepth, endDepth, startWidth, endWidth, falloffFactor);
        }
    }

    private static void CarveRiverSegment(float[,] accumulatedDepthMap, Point start, Point end, float startDepth,
        float endDepth, float startWidth, float endWidth, float falloffFactor)
    {
        // Calculate the line equation between start and end points
        Vector2 direction = (end - start).ToVector();
        int numSteps = (int)direction.Length();
        direction /= numSteps;

        Vector2 perpendicular = new Vector2(-direction.Y, direction.X);

        // Iterate over the points along the line between start and end points
        for (int i = 0; i <= numSteps; i++)
        {
            float t = (float)i / numSteps;

            Point position = start + Point.FromVector(direction * i);
            float currentDepth = Lerp(startDepth, endDepth, t);
            float currentWidth = Lerp(startWidth, endWidth, t);

            int widthInt = (int)Math.Ceiling(currentWidth);

            // Iterate over the points in the width direction
            for (int xOffset = -widthInt; xOffset <= widthInt; xOffset++)
            {
                for (int yOffset = -widthInt; yOffset <= widthInt; yOffset++)
                {
                    Point offsetPosition = position + new Point(xOffset, yOffset);

                    if (!IsInBounds(offsetPosition, accumulatedDepthMap)) continue;

                    float distance = new Vector2(xOffset, yOffset).Length();

                    if (distance <= currentWidth)
                    {
                        float scaledDistance = distance / currentWidth;
                        float depthFalloff = (float)(1 - (3 * Math.Pow(scaledDistance, 2) - 2 * Math.Pow(scaledDistance, 3)));
                        
                        float adjustedDepth = currentDepth * depthFalloff;
                        accumulatedDepthMap[offsetPosition.X, offsetPosition.Y] =
                            Math.Max(accumulatedDepthMap[offsetPosition.X, offsetPosition.Y], adjustedDepth);
                    }
                }
            }
        }
    }

    private static bool IsInBounds(Point position, float[,] map)
    {
        int width = map.GetLength(0);
        int height = map.GetLength(1);

        return position.X >= 0 && position.X < width && position.Y >= 0 && position.Y < height;
    }

    private static void ApplyAccumulatedDepthMap(float[,] heightMap, float[,] accumulatedDepthMap)
    {
        for (int x = 0; x < heightMap.GetLength(0); x++)
        {
            for (int y = 0; y < heightMap.GetLength(1); y++)
            {
                heightMap[x, y] = Math.Max(heightMap[x, y] - accumulatedDepthMap[x, y], 0);
            }
        }
    }

    private static float Lerp(float start, float end, float t)
    {
        return start + (end - start) * t;
    }
}