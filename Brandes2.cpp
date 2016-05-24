//Inspired from A Faster Algorithm for Betweenness Centrality paper.

//Takes O(nm) and O(nm + n^2log n) time on unweighted and weighted networks.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <list>

#include <cstring>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace std;

double clkbegin, clkend, t;

double rtclock(void)
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

// A simple neighbor struct, consisting of the target neighbor and the edge weight.
struct neighbor {
    int target;
    int weight;
    neighbor(int mTarget, int mWeight) : target(mTarget), weight(mWeight) {
    }
};

// The new adjacency list type.
typedef vector<vector<neighbor> > adjacency_list;

// Dijkstra's algorithm is used to calculate all the single source shortest paths in a weighted graph and the source's closeness.
float dijkstra_SSSP(int src, int n, stack<int> &visitStack, vector<int> &sigma, list<int> *pred, adjacency_list &adjList) {
    // Closeness counter.
    float closeness = 0;

    // Vector that holds the distances from the source.
    vector<int> dist;
    dist.resize(n, numeric_limits<int>::max());
    dist[src] = 0;
    
    // Queue used for the Dijkstra's algorithm.
    set<pair<int, int> > visitQueue;
    visitQueue.insert(make_pair(dist[src], src));

    // While there are still elements in the queue.
    while (!visitQueue.empty()) {
        // Pop the first.
        set<pair<int, int> >::iterator vPair = visitQueue.begin();
        int srcDist = vPair->first;
        int v = vPair->second;
        visitQueue.erase(vPair);
        visitStack.push(v);

        // Closeness part aggregation.
        closeness += srcDist;

        // Check the neighbors w of v.
        for (vector<neighbor>::iterator it = adjList[v].begin(); it != adjList[v].end(); it++) {
            int w = it->target;
            int newDist = srcDist + it->weight;
            // Node w found for the first time or the new distance is shorter?
            if (newDist < dist[w]) {
                visitQueue.erase(make_pair(dist[w], w));
                visitQueue.insert(make_pair(newDist, w));
                dist[w] = newDist;
                pred[w].clear();
                sigma[w] = 0;
            }
            // Is the shortest path to w via u?
            if (newDist == dist[w]) {
                pred[w].push_back(v);
                sigma[w] += sigma[v];
            }
        }
    }
    // Closeness part inversion.
    if (closeness!=0) {
        return 1.0 / closeness;
    } else {
        return 0;
    }
}

// BFS algorithm is used to calculate all the single source shortest paths in a non weighted graph and the source's closeness.
float bfs_SSSP(int src, int n, stack<int> &visitStack, vector<int> &sigma, list<int> *pred, adjacency_list &adjList) {
    // Closeness counter.
    float closeness = 0;
    
    // Vector that holds the distances from the source.
    vector<int> dist;
    dist.resize(n, -1);
    dist[src] = 0;

    // Queue used for the Bfs algorithm.
    queue<int> visitQueue;
    visitQueue.push(src);
    //cout<<"debug 30"<<endl;
    // While there are still elements in the queue.
    while (!visitQueue.empty()) {
        // Pop the first.
        int v = visitQueue.front();
        visitQueue.pop();
        visitStack.push(v);
        // Closeness part aggregation.
        closeness += dist[v];
        //cout<<"debug 31.1.0: "<<v<<endl;
        // Check the neighbors w of v.
        for (vector<neighbor>::iterator it = adjList[v].begin(); it != adjList[v].end(); it++) {
            int w = it->target;
            if(w < n && w > 0){
                // Node w found for the first time?
                if (dist[w] < 0) {
                    visitQueue.push(w);
                    dist[w] = dist[v] + 1;
                }
                //if(v==1)cout<<"debug 31.1.2: "<<v<<" "<<w<<endl;
                // Is the shortest path to w via u?
                if (dist[w] == dist[v] + 1) {
                    pred[w].push_back(v);
                    sigma[w] += sigma[v];
                }
            }
        }

    }    
    // Closeness part inversion.
    if (closeness!=0) {
        return 1.0 / closeness;
    } else {
        return 0;
    }
}

// Given two node indices, this function returns a string representation with
// the smallest index first, a dash, and the second index after.
string getEdgeTag(int n1, int n2) {
    ostringstream os;
    if (n1 <= n2) {
        os << n1 << "-" << n2;
    } else {
        os << n2 << "-" << n1;
    }
    return os.str();
}

// Prompts the user to type an input filename and returns the file pointer.
FILE * readPrompt() {
    FILE * fp;
    char str[100];
    
    // Prompt filename input.
    cout << "Enter filename: ";
    cin >> str;
    
    // Append .net if the user did not provide it.
    if (strstr(str, ".net") == NULL) {
        strcat(str, ".net");
    }
    
    // Open the file and return fp, if it exists, otherwise exit with error.
    fp = fopen(str, "r");
    if (fp == NULL) {
        cout << "File '" << str << "' not found.";
        exit(EXIT_FAILURE);
    } else {
        return fp;
    }
}

// Prints input file statistics just after input has finished.
void printInputStats(bool isWeigthed, int n, int e) {
    ofstream out;
    out.open ("out_graph_stats.txt");
    cout << "\n==================="
            << "\nINPUT GRAPH STATS"
            << "\n>Weighted: " << boolalpha << bool(isWeigthed)
            << "\n>#ofNodes: " << n
            << "\n>#ofEdges: " << e
            << "\n===================\n\n";
    out << "Weighted: " << boolalpha << bool(isWeigthed)
            << "\n>#ofNodes: " << n
            << "\n>#ofEdges: " << e;
    out.close();
}

// Reads an input file and fills up the adjacency list as well as the edges.
void readGraph(int &n, bool &isWeigthed, adjacency_list &adjList) {

    int e = 0; // Total number of edges (for statistics).
    isWeigthed = false;

    char * line = NULL;
    size_t len = 0;
    FILE * fp = readPrompt();

    // Find n, the total number of nodes.
    if (getline(&line, &len, fp) != -1) {
                strtok(line, " ");
            n = atoi(strtok(NULL, " "));

    }
    cout<<"read number of nodes :"<<n<<endl;
    // Reserve n space for adjacency list. If it fails, n was not parsed.
    if (n) {
        adjList.reserve(n);
    } else {
        cout << "Malformed input. Number of nodes undefined.";
        exit(EXIT_FAILURE);
    }
    cout<<"debug 1\n";
    
    // Read the nodes and the edges, one by one, and fill up adjList and edgeBetweenness.
    int start, end, weight;
    while (getline(&line, &len, fp) != -1) {
        e += 1;
        start = atoi(strtok(line, " ")) -1;
        end = atoi(strtok(NULL, " ")) -1;
        weight = 1;//atoi(strtok(NULL, " "));
        // Check if the graph is weighted. If w<=0, the input is malformed
        if (weight > 1) {
            isWeigthed = true;
        } else if(weight<=0) {
            cout << "Malformed input. Edge w weight=0.";
            exit(EXIT_FAILURE);
        }

        adjList[start].push_back(neighbor(end, weight));
        adjList[end].push_back(neighbor(start, weight));
    }
    cout<<"debug 10\n";
    if (line) {
        free(line);
    }
    
    // Print statistics after reading.
    printInputStats(isWeigthed, n, e);
}

// Clears the variables or re-initializes to 0, so that they are ready for the next loop.
void resetVariables(int src, int n, list<int> *pred, vector<int> &sigma, vector<float> &delta) {
    for (int i = 0; i < n; i++) {
        pred[i].clear();
    }

    sigma.clear();
    sigma.resize(n, 0);
    sigma[src] = 1;

    delta.clear();
    delta.resize(n, 0);
}

// Prints Closeness Centrality.
void printCloseness( int n, vector<float> closeness, bool normalize) {
    float nrml = 1;
    if (normalize) {
        nrml = 1.0/(n - 1);
    } 
    ofstream out;
    out.open ("out_closeness.txt");
    cout << "> Closeness Centrality" << endl;    
    for (int i = 0; i < n; i++) {
        //cout << "Node " << i << ": " << closeness[i] / nrml << endl;
        out << "Node " << i << ": " << closeness[i] / nrml << endl;
    }
    out.close();
}

// Prints Node Betweenness Centrality.
void printNodeBetweenness( int n, vector<float> nodeBetweenness, bool normalize) {
    float nrml = 1;
    if (normalize) {
        nrml = (n - 1)*(n - 2);
    }
    ofstream out;
    out.open ("out_node_betweenness.txt");
    cout << endl << "> Node Betweenness Centrality" << endl;
    for (int i = 0; i < n; i++) {
        //cout << "Node " << i << ": " << nodeBetweenness[i] / nrml << endl;
        out << "Node " << i << ": " << nodeBetweenness[i] / nrml << endl;
    }
    out.close();
}

int main(void) {
    int n; // Number of nodes
    bool isWeigthed; // Weighted graph check.
    adjacency_list adjList; // Adjacency list.

    // Centrality measures.
    readGraph(n, isWeigthed, adjList);
    
    vector<float> nodeBetweenness;
    nodeBetweenness.resize(n, 0);
    vector<float> closeness;
    closeness.resize(n, 0);

    // Input is read, and values are set to all the arguments.
    
    list<int> pred[n]; // List of predecessors of node v.
    vector<int> sigma;
    vector<float> delta;
    stack<int> visitStack; // Stack that holds the inverse order of visited nodes.
    clkbegin = rtclock();

    // For each node of the graph.
    
    for (int src = 0; src < n; src++) { 
        // Prepare the variables for the next loop.
        resetVariables(src, n, pred, sigma, delta);
        closeness[src] = bfs_SSSP(src, n, visitStack, sigma, pred, adjList);
        // Get the inverse order of visited nodes.
        while (!visitStack.empty()) {
            int w = visitStack.top();
            visitStack.pop();
            
            // For each predecessors of node w, do the math!
            for (list<int>::iterator it = pred[w].begin(); it != pred[w].end(); it++) {
                int v = *it;
                float c = ((float) sigma[v] / (float) sigma[w]) * (1.0 + delta[w]);

                delta[v] += c;
            }
            // Node betweenness aggregation part.
            if (w != src) {
                nodeBetweenness[w] += delta[w];
            }
        }

    }
    clkend = rtclock();
    t = clkend-clkbegin;
    cout << "\n" ;
    cout << "Time Taken : " << t;
    cout <<" Approx GFlops : " << 1.0*n*n/t/1e9 << "\n";
    cout << "\n";
    // Printing output.
    printCloseness(n, closeness, true);
    printNodeBetweenness(n, nodeBetweenness, true);
    
    return 0;
}
