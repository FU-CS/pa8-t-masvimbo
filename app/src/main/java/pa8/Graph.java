package pa8;

/** 
 * Graph interface for a directed, connected graph with vertices numbered from 0 to n-1.
 * For weighted graphs, the default weight of each edge is assumed to be 1.
 */
interface Graph {

    /**
     * Add an edge between two vertices.
     * @param v vertex 1 (0-indexed)
     * @param w vertex 2 (0-indexed)
     * 
     * When called for weighted graphs, the default weight of the edge is assumed to be 1.
     */
    public void addEdge(int v, int w);

    /**
     *  Add a weighted edge between two vertices.
     *  @param v vertex 1 (0-indexed)
     *  @param w vertex 2 (0-indexed)
     *  @param weight the weight of the edge
     * 
     * When called for unweighted graphs, it should be equivalent to addEdge(v, w).
     * 
     */
    public void addWeightedEdge(int v, int w, int weight);

    /** 
     * Perform a Breadth-First Search (BFS) starting from the specified vertex.
     * @param start the starting vertex
     * @return a String representing the order of vertices visited, e.g., "0 1 2".
     */ 
    public String bfs(int start);

    /** 
     * Perform a Depth-First Search (DFS) starting from the specified vertex.
     * @param start the starting vertex
     * @return a String representing the order of vertices visited, e.g., "0 2 1".
     */
    public String dfs(int start);

    /**
     * Determine if there is a cycle in the graph.
     * @return true if the graph contains a cycle, false otherwise
     */ 
    public boolean hasCycle();

    /**
     * Find the shortest path between two vertices.
     * Assumes the graph is unweighted (or weights are all 1).
     * 
     * @param v source vertex
     * @param w destination vertex
     * @return a String representing the shortest path, e.g., "0 -> 1 -> 2".
     *         If no path exists, return an empty string or a message indicating no path.
     */ 
    public String shortestPath(int v, int w);

}