using HydrologyMaps;
using Voronoi2;
using Point = HydrologyMaps.Point;

namespace RancHydrologyMaps;

public class HydrologyMapGen
{
    private readonly HydrologyParameters _parameters;
    Random _random = new Random();

    // Not currently used but can be useful for debugging algorithm. These edges will be drawn on the final map
    public static List<IEdge> DebugGraphEdges1 = new();
    public static List<IEdge> DebugGraphEdges2 = new();

    public HydrologyMapGen(HydrologyParameters parameters)
    {
        this._parameters = parameters;
    }
    
    public HydrologyMap GenerateIsland(int width, int seed, int stepIndex)
    {
        if (seed == -1)
        {
            seed = _random.Next();
            Console.WriteLine("Seed: " + seed);
        }

        _random = new Random(seed);

        float[,] heightmap = new float[width, width];

        // Gets a simple shape of an island where each pixel is 0 for ocean and 1 for land. The shape is assumed to have no lakes and the ocean no extra islands.
        HydrologyTerrainFormer.GetBaseHeightmapMask(heightmap, width, _random, _parameters);


        // Maintain a structure of points along the border of the island, used to determine how far from coast points/pixels are
        KdTree<Point> borderGraphNodeKdTree = new(
            point => point.X,
            point => point.Y
        );

        // Prepare to generate river graph
        HydrologyRiverGraphGenerator hydrologyRiverGraphGenerator = new(_parameters, _random);
        
        // Trace the shape of the island shape and find potential spots for river mouths, also fill in our borderGraphNodeKdTree
        (List<DirectedNode> riverMouthCandidates, List<GraphNode> extraVoronoiPoints) =
            hydrologyRiverGraphGenerator.GetRiverMouthCandidates(heightmap, _parameters.SpaceBetweenRiverMouthCandidates, borderGraphNodeKdTree);
        
        // Now generate a more realistic noise-based island heightmap, we use borderGraphNodeKdTree to slope terrain at border for beach-like structure
        HydrologyTerrainFormer.InterpolateHeightMapSimple(heightmap, borderGraphNodeKdTree, width, 10,
            _parameters.MinLandHeight);
        
        if (riverMouthCandidates.Count == 0)
        {
            Console.WriteLine("Error: Failed to find river mouths");
            return new HydrologyMap(heightmap, riverMouthCandidates, new List<RiverEdge>(), new GridCell[0, 0],
                new List<GraphEdge>(),
                new List<IGraphNode>());
        }
        
    
        // Take the river mouths and expand them inland into a graph
        (List<DirectedNode> allRiverNodes, List<RiverEdge> riverEdges) =
            hydrologyRiverGraphGenerator.ExpandRiverGrid(heightmap, riverMouthCandidates, stepIndex);
        
        // Make a new kd-tree to contain all river nodes in the river graph
        KdTree<GraphNode> riverNodesKdTree=  new(
            point => point.X,
            point => point.Y
        );

        allRiverNodes.ForEach(n => riverNodesKdTree.Insert(n));

        List<GraphEdge>? voronoiEdges = null;
        List<IGraphNode>? allPointsForVoronoi = null;
        if (_parameters.CalculateVoronoi)
        {
            allPointsForVoronoi = new List<IGraphNode>();
            allRiverNodes.ForEach(x => allPointsForVoronoi.Add(x));
            extraVoronoiPoints.ForEach(x => allPointsForVoronoi.Add(x));

            voronoiEdges = IslandVoronoi.GenerateVoronoiEdges(allPointsForVoronoi);
            voronoiEdges = IslandVoronoi.GenerateExtraEdges(riverMouthCandidates, voronoiEdges, heightmap);
        }
        
        allRiverNodes.RemoveAll(x => x.Type == NodeType.Border); // Not quite sure if this is still needed...
        
        allRiverNodes.ForEach(x => riverNodesKdTree.Insert(x));

        // A data structure that was meant to hold metadata for each pixel, not currently used but might be useful in the future.
        GridCell[,] gridCells = new GridCell[width, width];
        
        // Now we want to calculate the flow of water at each river node
        hydrologyRiverGraphGenerator.FindFlowRateForAllNodes(heightmap, riverNodesKdTree, width, allRiverNodes);
        
        // At the moment we don't require/assume that each node is at waterlevel and we might have "rivers" that run up in mountains that may or may not contain water but will at least carve paths
        hydrologyRiverGraphGenerator.GenerateRiverNodeHeights(riverMouthCandidates);
        
        // Carve the river graph onto the heightmap
        HydrologyTerrainFormer.CarveRivers(heightmap, riverMouthCandidates, _parameters.BaseRiverWidth, _parameters.RiverCarveDepth,
            _parameters.RiverCarveBankSlope, _parameters.RiverCarveFlowWidthInfluence);

        return new HydrologyMap(heightmap, allRiverNodes, riverEdges, gridCells, voronoiEdges, allPointsForVoronoi);
    }
}