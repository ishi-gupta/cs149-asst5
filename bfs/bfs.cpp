#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <vector>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
    
    // Use OpenMP to parallelize the clearing of the `present` array
    #pragma omp parallel for
    for (int i = 0; i < list->max_vertices; i++) {
        list->present[i] = false;
    }
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    list->present = (bool*)malloc(sizeof(bool) * list->max_vertices);
    vertex_set_clear(list);
}

#define THREADS 8
// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    //GO THROUGH EVERY NODE IN OUR CURRENT FRONTIER
    //#pragma omp parallel for 
    // #define THREADS 4
    #pragma omp parallel for schedule(dynamic, 100)
    //for each node id in the current frontier 
    for (int i=0; i<frontier->count; i++) {
        //get the node id for each node in the frontier
        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        //if node is the last node in the graph (node == g.num_nodes - 1). Getting indices for outgoin, dense array stuff!
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges // The last edge in the entire graph
                           : g->outgoing_starts[node + 1]; // Start of the next node's edges

        // attempt to add all neighbors to the new frontier
        // #pragma omp parallel for 
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            //IF NOT VISITED
            //SOMETHING ABOUT THIS VARIABLE DECLARATION IS SUS 
            int dist = -1; 
            #pragma omp atomic read
            dist = distances[outgoing];
            //test and test and set
            if (dist == NOT_VISITED_MARKER) { //reflect the lowest possible distance 
                if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1)) {  
                    //if we are visiting this neigbout from the first time, do the following
                    int index = 0;
                    #pragma omp atomic capture
                    index = new_frontier->count++;
                    new_frontier->vertices[index] = outgoing;
                }
            }
        }
    }
}

void print_graph(struct graph* g) {
    printf("Graph has %d nodes and %d edges.\n", g->num_nodes, g->num_edges);

    printf("Outgoing edges:\n");
    for (int node = 0; node < g->num_nodes; node++) {
        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1) 
                       ? g->num_edges 
                       : g->outgoing_starts[node + 1];

        printf("Node %d:", node);
        for (int edge = start_edge; edge < end_edge; edge++) {
            printf(" %d", g->outgoing_edges[edge]);
        }
        printf("\n");
    }

    printf("Incoming edges:\n");
    for (int node = 0; node < g->num_nodes; node++) {
        int start_edge = g->incoming_starts[node];
        int end_edge = (node == g->num_nodes - 1) 
                       ? g->num_edges 
                       : g->incoming_starts[node + 1];

        printf("Node %d:", node);
        for (int edge = start_edge; edge < end_edge; edge++) {
            printf(" %d", g->incoming_edges[edge]);
        }
        printf("\n");
    }
}


// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    printf("Running top-down BFS\n");
    //print_graph(graph);

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

int stepcount = 0;
void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances,
    std::vector<int>& remaining_nodes) // Pass the set of unvisited nodes)
{
    
    // printf("doing step %d\n", stepcount);
    stepcount++;

    double start_time = CycleTimer::currentSeconds();
    //j is the node we are checking to see if is in the fronteir
    int counter = 0;
    //TODO figure out what ynamic 100 represents 
    #pragma omp parallel for schedule(dynamic, 100) reduction(+:counter)//unifies thread local counters at the end 
    for (size_t i = 0; i < remaining_nodes.size(); i++){
        //printf("hi im checking node %d\n", node);
        int node = remaining_nodes[i];
        //TODO: no need to check if its not visited, we are only checking the remaining nodes
        if (distances[node] == NOT_VISITED_MARKER) {
            //check if it has an incoming edge from a node in the fronteir 
            // double inner_start_time = CycleTimer::currentSeconds();
            // for (int i=0; i<frontier->count; i++) {
                //collecting all the nodes incoming edges 
                int start_edge = g->incoming_starts[node];
                int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->incoming_starts[node + 1];
                //going through all the nodes who are incoming to me 
                double inner_start_time2 = CycleTimer::currentSeconds();
                for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                    int incoming = g->incoming_edges[neighbor];
                    if (frontier->present[incoming]) {
                        //TODO: thibk about wethe ur need compare and swap here 
                        if (distances[node] == NOT_VISITED_MARKER) {
                             distances[incoming] + 1;
                             new_frontier->present[node] = true;
                             counter++;
                            //think about why you need it 
                            //how do i restrucutre so I don't need this anymore 
                        }
                        // if (__sync_bool_compare_and_swap(&distances[node], NOT_VISITED_MARKER, distances[incoming] + 1)) { //my distance is the incoming edges already collected distance + 1  
                        //     //int index = 0;
                        //     new_frontier->present[node] = true;
                        //     // frontier->vertices[counter] = node; 
                        //     counter++; 
                        //     break;
                        //     // #pragma omp atomic capture
                        //     // counter = new_frontier->count++; //add this node to the new fonteir
                        //     // new_frontier->vertices[index] = node;
                        // }
                    }
                // }
                double inner_end_time = CycleTimer::currentSeconds();
            }
            
        }
    }
    double end_time = CycleTimer::currentSeconds();
    // printf("Outerloop time for step %d:  %.4f sec. %f percent of total time. \n", stepcount, end_time - start_time, (end_time - start_time)/end_time - start_time *100);
    // printf("Middleloop time for step %d:  %.4f sec %f percent of total time. \n", stepcount, inner_end_time - inner_start_time, (inner_end_time - inner_start_time)/(end_time - start_time) *100);
    // printf("Innerloop time for step %d:  %.4f sec %f percent of total time. \n", stepcount, inner_end_time - inner_start_time2,(inner_end_time - inner_start_time2)/(end_time - start_time) *100);


    //counting how many nodes are there in the next fronteir
    new_frontier->count = counter;

    // Update remaining nodes
    std::vector<int> updated_remaining_nodes;
    for (int node : remaining_nodes) {
        if (distances[node] == NOT_VISITED_MARKER) {
            updated_remaining_nodes.push_back(node);
        }
    }
    remaining_nodes = std::move(updated_remaining_nodes);


}


void bfs_bottom_up(Graph graph, solution* sol)
{

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    std::vector<int> remaining_nodes;
    for (int i = 0; i < graph->num_nodes; i++) {
        remaining_nodes.push_back(i);
    }
    // setup frontier with the root node
    //frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    frontier->present[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;
    frontier->count = 1;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step(graph, frontier, new_frontier, sol->distances, remaining_nodes);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
