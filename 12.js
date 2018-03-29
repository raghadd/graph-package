// CPCS 324 Algorithms & Data Structures 2
// Graph data structure starter - Final Project (NEW)
// 2017, Dr. Muhammad Al-Hashimi
// -----------------------------------------------------------------------
// simple graph object with linked-list edge implementation and minimal fields
// extra vertex and edge property fields to be added later as needed
//
var v = [],
    e = [];
// -----------------------------------------------------------------------
function main_graph()
{
    // create Heap object
    var heap = new Heap();
    heap.makeMin(false);
    heap.Hinsert(2, "a");
    heap.Hinsert(9, "b");
    heap.Hinsert(7, "c");
    heap.Hinsert(6, "d");
    heap.Hinsert(5, "e");
    heap.Hinsert(8, "f");
    // Print the heap  
    document.write(heap.show());

    heap.Hinsert(10, "g");
    document.write(heap.show());

    heap.Hinsert(15, "h");
    document.write(heap.show());

    /*heap.Hinsert(10,"i");
    document.write(heap.show());

    heap.Hinsert(4,"j");
    document.write(heap.show());

    heap.Hinsert(4,"k");
    document.write(heap.show());
    var deleted;
    deleted=heap.deleteRoot();
    document.write(deleted[0] ,", ", deleted[1]);
    // Print the heap 
    document.write(heap.show());
    
    deleted=heap.deleteRoot();
    document.write(deleted[0] ,", ", deleted[1]);
    // Print the heap 
    document.write(heap.show());

    deleted=heap.deleteRoot();
    //print the deleted node
    document.write(deleted[0] ,", ", deleted[1]);
    // Print the heap    
    document.write(heap.show());
    */

    // make graph of Exercise 9.2: 1b  
    var g = new Graph();
    g.label = "Exercise 9.2: 1b (Levitin, 3rd edition)";
    g.digraph = false;
    g.weighted = true;
    g.readGraph(v, e);
    g.printGraph();

    //preform prim with heap based priority queue

    MS = g.prim3();
    document.write("<p>", "MST by Prim2 (heap PQ)<br>");
    document.write("(-,0)", ", ");
    for (var i = 0; i < MS.length; i++)
        document.write("  (", MS[i].parent, ", ", MS[i].vtree, ") ");

    var g2 = new Graph();
    g2.digraph = false;
    g2.weighted = true;
    g2.readGraph(v, e);

    //preform first prim  
    g2.prim();
    chngtonumber(g2);
    document.write("<p>", "MST by first prim<br>");
    document.write("(-,0),");
    for (var i = 1; i < g2.Prim_Edges.length; i++)
    {
        // print edges info 
        document.write("  (", g2.Prim_Edges[i].v.label, ", ",
            g2.Prim_Edges[i].u.label, ")");
    }
}

// -----------------------------------------------------------------------
// similar to starters 11 (REMOVE network object)

function Vertex(v)
{
    // property fields
    this.label = v.label;
    this.visit = false;
    this.adjacent = new List();
    // -------------------- 
    // member methods   
    this.adjacentByID = adjacentByIdImpl;
    this.vertexInfo = vertexInfoImpl;
    this.insertAdjacent = insertAdjacentImpl;
    this.incidentEdge = incidentEdgeImpl;

}

// -----------------------------------------------------------------------
// similar to starter 11
/**
 * @constructor
 * Edge method
 * @param {intger} vert_i 
 * @param {integer} weight 
 */
function Edge(vert_i, weight)
{
    // property fields  
    this.target_v = vert_i; //Id of edge target vertex
    this.weight = weight === undefined ? null : weight; //Edge weight or cost.  
    this.flow = 0; //flow for argument edge i,j

}

// -----------------------------------------------------------------------
// similar to starter 11
/**
@constructor
Create graph 
*/
function Graph()
{
    // property fields
    this.vert = [];
    this.nv = 0;
    this.ne = 0; // number of edges
    this.digraph = false; // true if digraph, false otherwise (default undirected)
    this.dfs_push = []; // DFS traversal order output array
    this.bfs_order = []; // BFS traversal order output array
    this.label = ""; // identification string to label grap
    this.connectedComp = 0; // number of connected comps set by DFS; 0 (default) for no info
    this.AdjMatrix = []; // graph adjacency matrix to be created on demand
    this.weighted = false //true if weighted, false otherwise (default unweighted)
    /**#@+ 
     * @description Property field -
     */

    /** Counter for topoSearch chooses between DFS/BFS methods
    (1: dfs ) increase by 1 before bfs call
    @default 1 */
    this.topoCounter = 1;
    /** Transitive closure matrix , Created on demand,
    set after applying dfsTC. Stores output of last dfsTC call
    @default [ ] */
    this.TCMatrix = [];
    /** Transitive closure matrix , Created on demand,
    set after applying warshallFloyd. Stores output of last warshallFloyd call.
    @default [ ] */
    this.warshallTC = [];
    /** Distance matrix of shortest paths, Created on demand,
    set after applying warshallFloyd. Stores output of last warshallFloyd call.
    @default [ ] */
    this.floydD = [];

    /*array of distance id, Stores output of length
    of shortest path from source to current vertex*/
    this.dv = [];
    //Array of vertices 
    this.VT = [];
    //array of parent vertex id,Stores output of next to last vertex in path.
    this.parent = [];

    /** List of weighted edges making MST tree of graph, Created on demand,
    set after  prim. Stores output of last prim call. 
    @config {string} label Expected: vertex label
    @config {string} label Expected: vertex label
    @param {integer} w Expected: edge weight*/
    this.Prim_Edges = new Array();

    /**#@- */

    // member methods
    this.make_graph = make_graphImpl;
    this.read_graph = better_input; // ... (complete next)
    this.addEdge = addEdgeImpl2;
    this.readGraph = better_input; // default input reader method   
    this.printGraph = printGraphImpl; // better printer function
    this.list_vert = list_vert;
    this.componentInfo = componentInfoImpl //Get connectivity info in printable form
    this.isConnected = isConnectedImpl; //Return connected status based on internal connectivity field
    this.makeAdjMatrix = makeAdjMatrixImpl2; // initialize the adjacency matrix
    this.topoSearch = topoSearchImpl; // perform a topological search    
    this.dfs = dfsImpl; // DFS a connected component
    this.bfs = bfsImpl; // BFS a connected component    

    // --------------------     
    // transitive closure package (requirements in line comments) 
    // student methods next; actual functions in student code sections

    /**#@+
     * @description <b>Member method</b> -
     */

    // transitive closure package

    /** test if two vertices (v_i, v_j) have path in digraph */
    this.hasPath = hasPathImpl;
    /** return distance of shortest path between v_i, v_j in weighted graph   */
    this.shortestPath = shortestPathImpl;
    /** test if directed  acyclic */
    this.isDAG = isDAGImpl;
    /** create TC matrix , and distance matrix if weighted  */
    this.warshallFloyd = warshallFloydImpl;
    /** create DFS-Based TC matrix */
    this.dfsTC = dfsTCImpl;
    // --------------------  
    this.prim = primImpl;
    //return shortest path
    this.prim2 = prim2Impl;
    this.prim3 = prim3Impl;
    //single source shortest path
    //this.Dijkstra = DijkstraImpl;

    /**#@- */
}

//-----------------------------------------------------------------------
// functions used by methods of Graph and subsidiary objects

function make_graphImpl(n, m, w)
{
    // parameter validations and checks: number of edges etc.
    var mmax = n * (n - 1);
    if (!this.digraph) mmax /= 2;
    if (m > mmax)
    {
        document.write("<p>ERROR: invalid number of edges for graph type</p>");
        return;
    }

    // create n vertex in v[] using id 0 to n-1 as label
    var v = [];
    for (var i = 0; i < n; i++)
        v[i] = {
            label: i.toString()
        };

    // if graph complete no need to generate random edges, just create mmax edges systematically


    // otherwise repreat create m distinct edges (graph loops not allowed)

    var e = [],
        wmin = 1,
        wmax = 50000,
        wsum = 0;

    var h = []; // quick-dirty n x n matrix to check previously generated edges, 
    // m-entry hash table would be more efficient
    for (i = 0; i < n; i++)
    {
        h[i] = [];
        h[i][i] = 0; // no graph loops; 0 = blocked pair of vertices
    }

    for (i = 0; i < m; i++)
    {
        // generate vertices u, v randomly
        do {
            var u_i = random(0, n - 1),
                v_i = random(0, n - 1);

        } while (h[u_i][v_i] != undefined);

        h[u_i][v_i] = 0;
        h[v_i][u_i] = 0; // update matrix: block u,v; block v,u also if undirected

        // if (u,v) is distinct insert in e[] (generate random weight if w true)
        // otherwise repeat generate another u,v pair

        e[i] = {
            u: u_i,
            v: v_i
        };
        if (w)
        {
            e[i].w = random(wmin, wmax);
            wsum += e[i].w;
        }
    }

    // call graph reader method and set label, graph type depends on value of digraph property
    this.read_graph(v, e);
    this.label = "Generated " + n + " vertices, " + m + " random " + (!this.digraph ? "un" : "") + "directed edges (" + Math.round(m / mmax * 100) + "%)" + (w ? ", ave weight = " + Math.round(wsum / m) : "");
}

function random(low, high)
{
    return Math.floor(Math.random() * (high - low + 1)) + low;
}

// -----------------------------------------------------------------------
// begin student code section (REMOVE network functions)
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// implementation 1, transitive closure package + First Prim 
// transitive closure package 
/**
check if there is a path between two adjacency vertices by thier IDs
{@link #warshallTC} 

@methodOf Graph#
@param {Integer} v_i Source vertex id
@param {Integer} v_j target vertex id
@returns (boolean) true if there exists a path between v_i and v_j
*/
function hasPathImpl(v_i, v_j)
{
    return this.warshallTC[v_i][v_j] == 1 ? true : false;
}

// --------------------
/**
return the shortest path between two vertices by thier IDs      
{@link #floydD} 
@methodOf Graph#
@param {Integer} v_i Starting vertex id
@param {Integer} v_j end vertex id
@returns (Integer) the shortest path between v_i and v_j
*/
function shortestPathImpl(v_i, v_j)
{
    return this.floydD[v_i][v_j];
}

// --------------------
/**
Check if the main diagonal elements values equal zero,
to find that no cycles exixts.
{@link #hasPath}

@methodOf Graph#
@returns (boolean) true if caller graph is directed Acyclic graph
*/
function isDAGImpl()
{
    //return true if all the TC main diagonal values are zero 
    for (var i = 0, j = 0; i < this.warshallTC.length && j < this.warshallTC.length; i++, j++)
        if (this.hasPath(i, j))
        {
            return false;
        }
    return true;
}

// --------------------
/**
implementation of warshalFloyd algorithm 
Warshall's algorithm impelments the transitive closure matrix  {@link #warshallTC}
Floyd's algorithm finds the shortest path then impelments the distance matrix {@link #floydD}
@methodOf Graph#    

*/
function warshallFloydImpl()
{
    //  ADJACENCY MATRIX 
    this.makeAdjMatrix();

    //Fill warshallTC and distance matrix floydD arrays with adjacent matrix
    for (var i = 0; i < this.AdjMatrix.length; i++)
    {
        //take a copy of each row 
        this.warshallTC[i] = this.AdjMatrix[i].slice();
        this.floydD[i] = this.AdjMatrix[i].slice();
        for (var j = 0; j < this.nv; j++)
        { // when the value is 0 and i != j put infinity
            if (this.AdjMatrix[i][j] == 0 && i != j)
            {
                this.floydD[i][j] = Infinity;
            }
        }
    }

    // warshall-Floyed actual algorithm
    for (var k = 0; k < this.floydD.length; k++)
    {
        for (var i = 0; i < this.floydD.length; i++)
        {
            for (var j = 0; j < this.floydD.length; j++)
            {
                this.floydD[i][j] = Math.min(this.floydD[i][j], (this.floydD[i][k] + this.floydD[k][j]));
                this.warshallTC[i][j] = this.warshallTC[i][j] || (this.warshallTC[i][k] && this.warshallTC[k][j]) ? 1 : 0;
            }
        }
    }

    //replace Infinity with 0; because no distance
    for (var i = 0; i < this.floydD.length; i++)
        for (var j = 0; j < this.floydD.length; j++)
            if (this.floydD[i][j] == Infinity)
                this.floydD[i][j] = 0;
}

// --------------------

/**
dfs-based transitive clouser by going through each vertex and mark the visited by dfs with one
{@link #dfsTC}

@methodOf Graph# 

*/
function dfsTCImpl()
{
    for (var i = 0; i < this.nv; i++)
    {
        // get vertex i
        var v = this.vert[i];

        // mark all vertices unvisited
        for (var p = 0; p < this.nv; p++)
        {
            this.vert[p].visit = false;
        }

        // create and initilize each corresponding row 
        this.TCMatrix[i] = [];
        for (var j = 0; j < this.nv; j++)
            this.TCMatrix[i][j] = 0;

        // perform DFS search for each adjacent to the vertex v by its ID
        var w = v.adjacentByID();
        for (var n = 0; n < w.length; n++)
            this.dfs(w[n]); //for each adjacent vertex call dfs()

        //traverse the vertices to check which is visited
        for (var k = 0; k < this.nv; k++)
        {
            // if visited set 1 in the corresponding TC matrix
            if (this.vert[k].visit)
            {
                this.TCMatrix[i][k] = 1;
            }
        }
    }
}

/**
get the minimum spanning tree of a graph
{@link #prim}

@methodOf Graph#
@returns (list) list of MST edges of a graph 
*/
function primImpl()
{
    var Tvert = []; // declare the tree vertecies list

    var n;
    // mark all vertices unvisited    
    for (var w = 0; w < this.nv; w++)
    {
        this.vert[w].visit = false;
    }

    // get first vertex of graph and insert it at first index of Tree vertecies
    Tvert[0] = this.vert[0];
    // mark it as visited
    Tvert[0].visit = true;
    // create a record and init it with infinity ( used to compare distances)
    var record;
    // for all the vertecies 
    for (var i = 1; i < this.nv; i++)
    {
        record = 10000000000;
        // find minimum weight edge
        for (var j = 0; j < Tvert.length; j++)
        {
            // get incident edge of the vertex
            var incedint_Edge = Tvert[j].incidentEdge();
            // process all vertecies adjacent to the vert in tree vertecies
            // store the edge distance
            for (var k = 0; k < incedint_Edge.length; k++)
            {
                if (!this.vert[incedint_Edge[k].adjVert_i].visit &&
                    incedint_Edge[k].edgeWeight < record)
                {
                    this.Prim_Edges[i] = {
                        v: Tvert[j],
                        u: this.vert[incedint_Edge[k].adjVert_i],
                        w: incedint_Edge[k].edgeWeight
                    };
                    // update record 
                    record = this.Prim_Edges[i].w;
                }
            }
        }
        //add the adjacent vertex with min edge to Tree vertecies
        n = this.Prim_Edges.length;
        Tvert[Tvert.length] = this.Prim_Edges[n - 1].u;
        // mark it as visited
        this.Prim_Edges[n - 1].u.visit = true;
    }
    return this.Prim_Edges;
}
// -----------------------------------------------------------------------
// implementation 2, greedy algorithms package (REMOVE Dijkstra)
/**
 * 
 * Function that finds the MST of a weighted graph by Prim's algorithm's second implementation 
 * using a priority queue.
 * 
 * @methodOf Graph#
 * @author Raghad Salem
 * 
 * @returns {object[]} Array of custom objects containing MST edges infromation. 
 * Fields of objects are as follows:
 * 
 * @returns {integer} <code>object.parent</code> - id of source vertex
 * @returns {integer} <code>object.vtree</code> - id of destination vertex
 * @returns {integer} <code>object.key</code> - weight of edge
 * 
 */

function prim2Impl()
{
    // Create Priority Queue
    var que = new PQueue();
    // empty arrays
    var Vt = [],
        Et = [],
        // var to store incident edges'
        inci_edges,
        // var used to store dequeued objects
        deqVert;
    // Mark all verticies unvisited
    for (var i = 0; i < this.nv; i++)
        this.vert[i].visit = false;

    Vt[0] = this.vert[0];
    // insert first vertex in tree and mark visited
    this.vert[0].visit = true;

    // Adding incident edges to priority queue (first time)
    inci_edges = Vt[0].incidentEdge();
    for (var j = 0; j < inci_edges.length; j++)
        que.insert(
        {
            parent: 0,
            vtree: inci_edges[j].adjVert_i
        }, inci_edges[j].edgeWeight);

    for (var k = 1; k < this.nv; k++)
    {
        // find first unvisited vert (w/ shorest path) 
        deqVert = que.deleteMin();
        while (this.vert[deqVert.dataitem.vtree].visit)
            deqVert = que.deleteMin();
        // Add vertex to visited tree (array) 
        Vt[Vt.length] = this.vert[deqVert.dataitem.vtree];
        this.vert[deqVert.dataitem.vtree].visit = true;
        // Add edges to MST
        Et[Et.length] = {
            parent: deqVert.dataitem.parent,
            vtree: deqVert.dataitem.vtree,
            w: deqVert.key
        };
        // Adding new incident edges to priority queue
        inci_edges = Vt[(Vt.length - 1)].incidentEdge();
        for (var l = 0; l < inci_edges.length; l++)
            if (inci_edges[l].edgeWeight != 0)
                que.insert(
                {
                    parent: deqVert.dataitem.vtree,
                    vtree: inci_edges[l].adjVert_i
                }, inci_edges[l].edgeWeight);

    }

    return Et;

}

function prim3Impl()
{
    // Create Priority Queue
    var que = new HPQueue();

    // empty arrays
    var Vt = [],
        Et = [],
        // var to store incident edges
        inci_edges,
        // var used to store dequeued objects
        deqVert;
    // Mark all verticies unvisited
    for (var i = 0; i < this.nv; i++)
        this.vert[i].visit = false;

    Vt[0] = this.vert[0];
    // insert first vertex in tree and mark visited
    this.vert[0].visit = true;

    // Adding incident edges to priority queue (first time)
    inci_edges = Vt[0].incidentEdge();
    for (var j = 0; j < inci_edges.length; j++)
        que.hpinsert(
        {
            parent: 0,
            vtree: inci_edges[j].adjVert_i
        }, inci_edges[j].edgeWeight);

    for (var k = 1; k < this.nv; k++)
    {
        // find first unvisited vert (w/ shorest path) 
        deqVert = que.hpdeleteMin();
        while (this.vert[deqVert[1].vtree].visit)
            deqVert = que.hpdeleteMin();
        // Add vertex to visited tree (array) 
        Vt[Vt.length] = this.vert[deqVert[1].vtree];
        this.vert[deqVert[1].vtree].visit = true;
        // Add edges to MST
        Et[Et.length] = {
            parent: deqVert[1].parent,
            vtree: deqVert[1].vtree,
            w: deqVert[0]
        };
        // Adding new incident edges to priority queue
        inci_edges = Vt[(Vt.length - 1)].incidentEdge();
        for (var l = 0; l < inci_edges.length; l++)
            if (inci_edges[l].edgeWeight != 0)
                que.hpinsert(
                {
                    parent: deqVert[1].vtree,
                    vtree: inci_edges[l].adjVert_i
                }, inci_edges[l].edgeWeight);

    }

    return Et;

}

// -----------------------------------------------------------------------
// additional functions NOT in published API

// -----------------------------------------------------------------------
// paste your Heap() object, followed by functions implementing its methods

function Heap()
{
    // h[0] not used, heap initially empty

    this.h = [null]; // heap of integer keys
    this.h_item = [null]; // corresponding heap of data-items (any object)
    this.size = 0; // 1 smaller than array (also index of last child)
    this.isMinHeap = true;


    // --------------------
    // min pq-required; many more heap processing methods could be added here
    // the 2 basic shape maintainig operations heapify and reheapify simplify
    // processing functions

    this.isEmpty = lsEmptyimpl; // return true if heap empty
    this.deleteRoot = deleteRootImpl; // return data-item in root
    this.Hinsert = InsertImpl; // insert data-item with key
    this.makeMin = makeMinImpl;

    this.heapify = HeapifyImpl; // make subtree heap; top-down heapify ("sink") used by .deleteRoot()
    this.reheapify = reheapifyimpl; // bottom-up reheapify ("swim") used by .insert()
    this.show = heapShow; // utility: return pretty formatted heap as string 

    // --------------------
    // student methods next; ; actual functions in student code section at end


}

// -----------------------------------------------------------------------
// functions used by Heap() object methods
//

function heapShow()
{
    var n = this.size;
    var m = Math.floor(n / 2); // last parent node

    var k = this.h.slice(1, n + 1),
        a = this.h_item.slice(1, n + 1);

    var out = "<h2>Heap (size=" + n + "):</h2><p>Keys: " + k + "<br>Data: " + a + "</p>";
    for (var i = 1; i <= m; i++)
    {
        out += "<p>" + i + ": <b>" + this.h[i] + "(" + this.h_item[i] + ")</b><ul>";
        if (2 * i <= n)
            out += "<li>" + this.h[2 * i] + "</li>";
        if (2 * i + 1 <= n)
            out += "<li>" + this.h[2 * i + 1] + "</li>";
        out += "</ul></p>";
    }

    return out;
}
// ------------------------
/**
 * 
 * Function that checks if heap is empty or not
 * 
 * @methodOf Heap#
 * @author Nura derya
 * 
 * @returns {Boolean} true or false
 */
function lsEmptyimpl()
{
    return (this.size == 0);
}
// ------------------------
/**
 * 
 * Function that deletes the root whtch is the largest by swapping the root with the leaf node returns it
 * use heapify to fix any imbalance in the heap recursively
 * 
 * @methodOf Heap#
 * @author Nura derya
 * @author Hatoon Mohammed
 * 
 * @returns {array} [0] the key and [1] the item of the root
 */
function deleteRootImpl()
{
    var root = [this.h[1], this.h_item[1]]; //root is at 1


    if (!this.isEmpty())
    {
        // swap between root and leaf element and decrease the size
        this.h_item[1] = this.h_item[this.size];
        this.h[1] = this.h[this.size--];

        //this.size--;
        this.heapify(1);
    }
    return root;
}
// ------------------------
/**
 * 
 * Function that inserts a node into the heap
 * and use reheapify to fix any imbalance recursively 
 * @methodOf Heap#
 * @author Nura derya
 * @author Hatoon Mohammed
 * 
 *
 */
function InsertImpl(key, data)
{
    if (this.isMinHeap)
    {
        key *= -1;
    }
    this.h[++this.size] = key;
    this.h_item[this.size] = data; //add it to the heap

    // n++;//increase heap size
    this.reheapify(); //fix the condition (min heap)
}
// ------------------------
/**
 * 
 * Function that fix any imbalance in the heap 
 *  by swapping between the parent and the largest of it's children
 * @methodOf Heap#
 * @author Hatoon Mohammed
 * @author Raghad Salem
 * 
 */
function HeapifyImpl(i)
{
    /*var n = this.size;
    while ((2 * i) <= n)
    {
        var j = 2 * i;
        if (j < n && (j < (j + 1))) j++;
        if (!(i < j)) break;

        var v = this.h[j];
        var u = this.h_item[j];

        this.h[j] = this.h[i];
        this.h_item[j] = this.h_item[i];

        this.h[i] = v;
        this.h_item[i] = u;

        i = j;
    }*/
    var largest = i;
    var l = 2 * i; //left
    var r = l + 1 //right
    //if left child is larger than root
    if (l < this.h.length && this.h[l] > this.h[largest])
    {
        largest = l;
    }
    //if right child is  larger than root
    if (r < this.h.length && this.h[r] > this.h[largest])
    {
        largest = r;
    }
    //if largest is not in root 
    if (largest != i)
    {
        var swap = this.h[i];
        var swap_item = this.h_item[i];
        this.h[i] = this.h[largest];
        this.h_item[i] = this.h_item[largest];
        this.h[largest] = swap;
        this.h_item[largest] = swap_item;
        //recursevvly heapfy the affected subtree
        this.heapify(largest);
    }
}

// ------------------------
/**
 * 
 * Function that fix any imbalance in the heap 
 *  by swappimg between the node and it's parent if it's larger than the parent
 * @methodOf Heap#
 * @author Aisha Abdullah
 * @author Raghad Salem
 */
function reheapifyimpl()
{
    /*var k = this.size;

    while (k > 1 && (this.h[Math.floor(k / 2)] < this.h[k]))
    {
        var v = this.h[Math.floor(k / 2)];
        var u = this.h_item[Math.floor(k / 2)];

        this.h[Math.floor(k / 2)] = this.h[k];
        this.h_item[Math.floor(k / 2)] = this.h_item[k];

        this.h[k] = v;
        this.h_item[k] = u;

        k = Math.floor(k / 2);

    }*/
    var n = this.size;
    // last parent node
    var l = Math.floor(n / 2);

    for (var i = l; i >= 1; i--)
    {
        var k = i;
        var v = this.h[k];
        var u = this.h_item[k];
        var heap = false;

        while (!heap && 2 * k <= n)
        {
            var j = 2 * k;

            if (j < n)
            { //there are two children
                if (this.h[j] < this.h[j + 1])
                {
                    j++;
                }

            }
            if (v >= this.h[j])
            {
                heap = true;
            }
            else
            {
                this.h[k] = this.h[j];
                this.h_item[k] = this.h_item[j];
                k = j;
            }
        } //End while loop
        this.h[k] = v;
        this.h_item[k] = u;

    } //end for loop
}

function makeMinImpl(b)
{
    this.isMinHeap = b;
}
// ...etc

// -----------------------------------------------------------------------
// similar to starter 8 and 11
// published docs section (ref. assignment page)
// for this section, strip line comments (leave outline)
// *NO* JSDOC comments in this section
// -----------------------------------------------------------------------

function better_input(v, e)
{
    // set number of vertices and edges fields
    this.nv = v.length;
    this.ne = e.length;
    // input vertices into internal vertex array
    var i;
    for (i = 0; i < this.nv; i++)
    {
        this.vert[i] = new Vertex(v[i]);
    }
    // input vertex pairs from edge list input array
    // remember to pass vertex ids to add_edge() 
    for (var i = 0; i < this.ne; i++)
    {
        this.addEdge(e[i].u, e[i].v, e[i].w);
    }
    // double edge count if graph undirected 
    if (!this.digraph)
    {
        this.ne = e.length * 2;
    }
    // test if graph weighted
    if (!(e[0].w === undefined))
    {
        this.weighted = true;
    }
}
// --------------------
/**
 * @methodOf Graph#  
 */
function list_vert()
{
    return "";
}

// --------------------
/**
*   @methodOf Graph# 

*/
function printGraphImpl()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "" : "UN", "WEIGHTED, ",
        this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");

    for (var i = 0; i < this.nv; i++)
    {
        var v = this.vert[i];
        document.write("VERTEX: ", i, v.vertexInfo(), "<br>");
    }
    document.write("<br>");
}

// --------------------
/**
 * @methodOf Graph#  
 */
function DFS()
{

}

// --------------------
/**
 * @methodOf Graph#  
 */
function dfsImpl(v_i)
{
    // process vertex
    var v = this.vert[v_i];
    v.visit = true;
    this.dfs_push[this.dfs_push.length] = v_i;
    // recursively traverse unvisited adjacent vertices 
    var w = v.adjacentByID();
    var i;
    for (i = 0; i < w.length; i++)
    {
        if (!this.vert[w[i]].visit)
        {
            this.dfs(w[i]);
        }
    }

}

// --------------------
/**
 * @methodOf Graph#  
 */
function BFS()
{

}
// --------------------
/**
 * @methodOf Graph#  
 */
function bfsImpl(v_i)
{
    // get vertex v by its id
    var v = this.vert[v_i];

    // process v 
    v.visit = true;
    this.bfs_order[this.bfs_order.length] = v_i;

    // initialize queue with v
    var Q = new Queue();
    Q.enqueue(v);

    // while queue not empty
    while (!Q.isEmpty())
    {
        // dequeue and process a vertex, u
        var u = Q.dequeue();

        // queue all unvisited vertices adjacent to u
        var w = u.adjacentByID();
        var i;
        for (i = 0; i < w.length; i++)
        {
            if (!this.vert[w[i]].visit)
            {
                this.vert[w[i]].visit = true;
                Q.enqueue(this.vert[w[i]]);
                this.bfs_order[this.bfs_order.length] = w[i];
            }
        }
    }
}

// --------------------
/**
 * @methodOf Graph#  
 */
function addEdgeImpl(u_i, v_i)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex 
    var u = this.vert[u_i];
    var v = this.vert[v_i];

    // insert (u,v), i.e., insert v (by id) in adjacency list of u
    u.adjacent.insert(v_i);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        v.adjacent.insert(u_i);
    }

}

// --------------------
/**
 * @methodOf Graph#  
 */
function addEdgeImpl2(u_i, v_i, w)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.vert[u_i];
    var v = this.vert[v_i];
    // insert (u,v), i.e., insert v in adjacency list of u
    u.insertAdjacent(v_i, w);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        v.insertAdjacent(u_i, w);
    }

}

// --------------------  
/**
 * @methodOf Graph#  
 */
function topoSearchImpl(fun)
{
    // mark all vertices unvisited
    var i;
    for (i = 0; i < this.nv; i++)
    {
        this.vert[i].visit = false;
    }

    // traverse unvisited connected component 

    for (i = 0; i < this.nv; i++)
    {
        if (!this.vert[i].visit)
        {
            this.topoCounter == 1 ? (fun++, this.dfs(i)) : this.bfs(i);
        }
    }
    return fun;
}

// --------------------
function makeAdjMatrixImpl()
{

}

// --------------------
/**
 * @methodOf Vertex#  
 */
function adjacentByIdImpl()
{
    var adjacent_ID = [];
    var edges_obj = this.adjacent.traverse();
    for (var i = 0; i < edges_obj.length; i++)
    {
        adjacent_ID[i] = edges_obj[i].target_v;
    }

    return adjacent_ID;

}

//----------------------
/**
 *  @methodOf Vertex#  
 */
function insertAdjacentImpl(v_i, weight)
{
    var edge = new Edge()
    edge.target_v = v_i;
    edge.weight = weight;
    this.adjacent.insert(edge);
}

/**
 * @methodOf Vertex#  
 */
function vertexInfoImpl()
{
    return " {" + this.label + "} - VISIT: " + this.visit + " - ADJACENCY: " + this.adjacentByID();
}

//-----------------------
/**
 * @methodOf Graph#  
 */
function componentInfoImpl()
{
    //return CONNECTED, if graph is connected
    if (this.isConnected())
        return "CONNECTED";
    //return No connectivity info if there is no connected components
    else if (this.connectedComp == 0)
        return "no connectivity info";
    //return DISCONNECTED if  it has more than 1 connected component
    else
        return "DISCONNECTED " + this.connectedComp;
}

// --------------------
/**
 * @methodOf Graph#  
 */
function isConnectedImpl()
{
    return this.connectedComp == 1 ? true : false;
}

// --------------------
/**
 * @methodOf Graph#
 */
function makeAdjMatrixImpl2()
{
    // initially create row elements and zero the adjacncy matrix
    for (var i = 0; i < this.nv; i++)
    {
        var v = this.vert[i];
        this.AdjMatrix[i] = [];
        for (var j = 0; j < this.nv; j++)
        {
            this.AdjMatrix[i][j] = 0;
        }
        // process vertices adjacent to v: get incident edges
        // for each node, use target id to set matrix value
        var incident_Edge = v.incidentEdge();
        for (var k = 0; k < incident_Edge.length; k++)
        {
            // set adjacent target vertex by id
            var incident_Target = incident_Edge[k].adjVert_i;
            // a text label="" 
            var incident_Label = incident_Edge[k].edgeLabel;
            // weight/cost if any
            var incident_Weight = incident_Edge[k].edgeWeight;
            //set the adjacent vertex with the weight if any
            this.AdjMatrix[i][incident_Target] = this.weighted ? incident_Weight : 1;
        }
    }
}
/**
@methodOf Vertex# 
*/
function incidentEdgeImpl()
{
    var edgesInfo = []; // create and declare object arrays
    var w = this.adjacent.traverse(); //get incedint edges
    //get info of each edge in adjacent list
    for (var i = 0; i < w.length; i++)
    {
        edgesInfo[i] = {
            adjVert_i: w[i].target_v,
            edgeLabel: "",
            edgeWeight: w[i].weight
        };
    }
    return edgesInfo;
}

function chngtonumber(g)
{
    for (var i = 1; i < g.Prim_Edges.length; i++)
    {
        if (g.Prim_Edges[i].v.label == 'a')
        {
            g.Prim_Edges[i].v.label = 0;

        }
        else if (g.Prim_Edges[i].v.label == "b")
        {
            g.Prim_Edges[i].v.label = 1;

        }
        else if (g.Prim_Edges[i].v.label == "c")
        {
            g.Prim_Edges[i].v.label = 2;

        }
        else if (g.Prim_Edges[i].v.label == "d")
        {
            g.Prim_Edges[i].v.label = 3;

        }
        else if (g.Prim_Edges[i].v.label == "e")
        {
            g.Prim_Edges[i].v.label = 4;

        }
        else if (g.Prim_Edges[i].v.label == "f")
        {
            g.Prim_Edges[i].v.label = 5;

        }
        else if (g.Prim_Edges[i].v.label == "g")
        {
            g.Prim_Edges[i].v.label = 6;

        }
        else if (g.Prim_Edges[i].v.label == "h")
        {
            g.Prim_Edges[i].v.label = 7;

        }
        else if (g.Prim_Edges[i].v.label == "i")
        {
            g.Prim_Edges[i].v.label = 8;

        }
        else if (g.Prim_Edges[i].v.label == "j")
        {
            g.Prim_Edges[i].v.label = 9;

        }
        else if (g.Prim_Edges[i].v.label == "k")
        {
            g.Prim_Edges[i].v.label = 10;

        }
        else if (g.Prim_Edges[i].v.label == "l")
        {
            g.Prim_Edges[i].v.label = 11;

        }
        else
        {

        }

    }
    for (var i = 1; i < g.Prim_Edges.length; i++)
    {
        if (g.Prim_Edges[i].u.label == 'a')
        {

            g.Prim_Edges[i].u.label = 0;
        }
        else if (g.Prim_Edges[i].u.label == "b")
        {

            g.Prim_Edges[i].u.label = 1;
        }
        else if (g.Prim_Edges[i].u.label == "c")
        {

            g.Prim_Edges[i].u.label = 2;
        }
        else if (g.Prim_Edges[i].u.label == "d")
        {

            g.Prim_Edges[i].u.label = 3;
        }
        else if (g.Prim_Edges[i].u.label == "e")
        {

            g.Prim_Edges[i].u.label = 4;
        }
        else if (g.Prim_Edges[i].u.label == "f")
        {

            g.Prim_Edges[i].u.label = 5;
        }
        else if (g.Prim_Edges[i].u.label == "g")
        {

            g.Prim_Edges[i].u.label = 6;
        }
        else if (g.Prim_Edges[i].u.label == "h")
        {

            g.Prim_Edges[i].u.label = 7;
        }
        else if (g.Prim_Edges[i].u.label == "i")
        {

            g.Prim_Edges[i].u.label = 8;
        }
        else if (g.Prim_Edges[i].u.label == "j")
        {

            g.Prim_Edges[i].u.label = 9;
        }
        else if (g.Prim_Edges[i].u.label == "k")
        {

            g.Prim_Edges[i].u.label = 10;
        }
        else if (g.Prim_Edges[i].u.label == "l")
        {

            g.Prim_Edges[i].u.label = 11;
        }
        else
        {

        }

    }
}