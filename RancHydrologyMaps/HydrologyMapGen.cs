using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using NetTopologySuite.Geometries;
using NetTopologySuite.Index.KdTree;
using RancHydrologyMaps;
using Voronoi2;

namespace HydrologyMaps;

public class HydrologyMapGen
{
    private readonly HydrologyParameters parameters;
    Random Random = new Random();

    public static List<IEdge> DebugGraphEdges1 = new();
    public static List<IEdge> DebugGraphEdges2 = new();

    public HydrologyMapGen(HydrologyParameters parameters)
    {
        this.parameters = parameters;
    }


    public HydrologyMap GenerateIsland(int width, int seed, int k)
    {
        if (seed == -1)
        {
            seed = Random.Next();
            Console.WriteLine("Seed: " + seed);
        }

        Random = new Random(seed);


        float[,] heightmap = new float[width, width];

        HydrologyTerrainFormer.FillBaseHeightmap(heightmap, width, Random, parameters);


        KdTree<Point> borderGraphNodeKdTree = new(
            point => point.X,
            point => point.Y
        );


        HydrologyRiverGraph hydrologyRiverGraph = new(parameters, Random);
        
        
        (List<DirectedNode> riverMouthCandidates, List<GraphNode> extraVoronoiPoints) =
            hydrologyRiverGraph.GetRiverMouthCandidates(heightmap, parameters.SpaceBetweenRiverMouthCandidates, borderGraphNodeKdTree);
        
        
        HydrologyTerrainFormer.InterpolateHeightMapSimple(heightmap, borderGraphNodeKdTree, width, 10,
            5f);
        
        if (riverMouthCandidates.Count == 0)
        {
            return new HydrologyMap(heightmap, riverMouthCandidates, new List<RiverEdge>(), new GridCell[0, 0],
                new List<GraphEdge>(),
                new List<IGraphNode>());
        }


        KdTree<GraphNode> riverNodesKdTree=  new(
            point => point.X,
            point => point.Y
        );
        (List<DirectedNode> allRiverNodes, List<RiverEdge> riverEdges) =
            hydrologyRiverGraph.ExpandRiverGrid(heightmap, riverMouthCandidates, k);


        
        allRiverNodes.ForEach(n => riverNodesKdTree.Insert(n));

        List<IGraphNode> allPointsForVoronoi = new List<IGraphNode>();
        allRiverNodes.ForEach(x => allPointsForVoronoi.Add(x));
        extraVoronoiPoints.ForEach(x => allPointsForVoronoi.Add(x));

        List<GraphEdge> voronoiEdges = IslandVoronoi.GenerateVoronoiEdges(allPointsForVoronoi);
        voronoiEdges = IslandVoronoi.GenerateExtraEdges(riverMouthCandidates, voronoiEdges, heightmap);
        //   voronoiEdges = voronoiEdges.Where(x => Distance((int)x.x1, (int)x.x2, (int)x.y1, (int)x.y2) < 51.0).ToList();

        allRiverNodes.RemoveAll(x => x.Type == NodeType.Border);

        GridCell[,] gridCells = new GridCell[width, width];


        allRiverNodes.ForEach(x => riverNodesKdTree.Insert(x));

        hydrologyRiverGraph.FindAreaForNodes(heightmap, riverNodesKdTree, allRiverNodes, width, width);
        // flowrate = 0.42 · A^0.69 formula from hydrology paper. Before this point FlowRate holds the number of pixels for which the closest node is this
        allRiverNodes.ForEach(n => n.FlowRate = 0.42 * Math.Pow(n.FlowRate, 0.69));
        riverMouthCandidates.ForEach(n => hydrologyRiverGraph.SetCumulativeFlowRate(n));

        hydrologyRiverGraph.GenerateRiverNodeHeights(riverMouthCandidates);


        HydrologyTerrainFormer.CarveRivers(heightmap, riverMouthCandidates /*.Skip(1).Take(1).ToList()*/, 0.5f, 0.35f,
            0.3f, 1f, 8f);

        return new HydrologyMap(heightmap, allRiverNodes, riverEdges, gridCells, voronoiEdges, allPointsForVoronoi);
    }
}