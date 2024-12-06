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

void vertex_set_clear_bottom(vertex_set* list) {
    list->count = 0;
    
    // Use OpenMP to parallelize the clearing of the `present` array
    #pragma omp parallel for
    for (int i = 0; i < list->max_vertices; i++) {
        list->present[i] = false;
    }
}

void vertex_set_clear_top(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init_bottom(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    list->present = (bool*)malloc(sizeof(bool) * list->max_vertices);
    vertex_set_clear_bottom(list);
}

void vertex_set_init_top(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    list->present = (bool*)malloc(sizeof(bool) * list->max_vertices);
    vertex_set_clear_top(list);
}

//takes a frontier in top down format and gets it ready for bottom up bfs
void vertex_set_top_to_bottom(vertex_set* list){
    // #pragma omp parallel
    // printf("list count: %d\n", list->count);
    // printf("vetices array:\n"); 
    // for (int i = 0; i < list->count; i++){
    //     printf("%d\n", list->vertices[i]);
    // }
    #pragma omp parallel for
    for (int i = 0; i < list->count; i++){
        list->present[list->vertices[i]] = true;
    }
}

// void print_frontier(vertex_set* list){
//     printf("frontier count: %d\n", list->count);
//     printf("frontier vertices: \n");
//     for (int i = 0; i < list->count; i++){
//         printf("%d\n", list->vertices[i]);
//     }
//     printf("frontier present: \n");
//     for (int i = 0; i < list->max_vertices; i++){
//         printf("%d\n", list->present[i]);
//     }
// }

//takes a frontier in bottom up format and gets it ready for top down bfs
void vertex_set_bottom_to_top(vertex_set* list){
    // #pragma omp parallel
    int pos = 0; 
    for (int i = 0; i < list->max_vertices; i++){
        if (list->present[i]){
            list->vertices[pos] = i;
            pos++;
        }
    }
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    // printf("frontier\n");
    // // print_frontier(frontier);
    // printf("new frontier\n");
    // print_frontier(new_frontier);

    //GO THROUGH EVERY NODE IN OUR CURRENT FRONTIER
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




// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    printf("Running top-down BFS\n");
    //print_graph(graph);

    vertex_set list1;
    vertex_set list2;
    vertex_set_init_top(&list1, graph->num_nodes);
    vertex_set_init_top(&list2, graph->num_nodes);

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

        vertex_set_clear_top(new_frontier);

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




void bottom_up_step_2(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    double start_time = CycleTimer::currentSeconds();
    
    int counter = 0;
    
    #pragma omp parallel for schedule(dynamic, 100) reduction(+:counter)
    for (size_t node = 0; node < g->num_nodes; node++){
        if (distances[node] == NOT_VISITED_MARKER) {
            //check if it has an incoming edge from a node in the fronteir 
            int start_edge = g->incoming_starts[node];
            int end_edge = (node == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[node + 1];
            //going through all the nodes who are incoming to me
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++) {
                int incoming = g->incoming_edges[neighbor];
                if (frontier->present[incoming]) {
                    //TODO: thibk about wethe ur need compare and swap here 
                    // if (__sync_bool_compare_and_swap(&distances[node], NOT_VISITED_MARKER, distances[incoming] + 1)) { //my distance is the incoming edges already collected distance + 1  
                    //     new_frontier->present[node] = true;
                    //     counter++; 
                    //     break;
                    // }
                        if (distances[node] == NOT_VISITED_MARKER) {
                             distances[node] = distances[incoming] + 1;
                             new_frontier->present[node] = true;
                             counter++;
                            //think about why you need it 
                            //how do i restrucutre so I don't need this anymore 
                        }
                }
            }
        }
    }
    
    //counting how many nodes are there in the next fronteir
    new_frontier->count = counter;
}


void bfs_bottom_up(Graph graph, solution* sol)
{

    vertex_set list1;
    vertex_set list2;
    vertex_set_init_bottom(&list1, graph->num_nodes);
    vertex_set_init_bottom(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    
    // setup frontier with the root node
    frontier->present[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;
    frontier->count = 1;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear_bottom(new_frontier);

        bottom_up_step_2(graph, frontier, new_frontier, sol->distances);

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
    printf("starting bfs hybrid\n");
    //run top down step 
    vertex_set list1;
    vertex_set list2;
    vertex_set_init_bottom(&list1, graph->num_nodes);
    vertex_set_init_bottom(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node to start with top down
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    // printf("count: %d\n", frontier->count);

    //setup frontier with the root node to start with bottom up
    frontier->present[ROOT_NODE_ID] = true;


    // // swap pointers
    // vertex_set* tmp = frontier;
    // frontier = new_frontier;
    // new_frontier = tmp;

    int total_iterations = 0; 
    int top_down_iterations = 0;    
    int flag = 1; 
    while (frontier->count != 0) {
        //printf("frontier count: %d\n", frontier->count);

        if (flag ==1){
            vertex_set_clear_top(new_frontier);
            top_down_step(graph, frontier, new_frontier, sol->distances); 

            //printf("ratio = %f\n", ((float)new_frontier->count/ graph->num_nodes));
            if (((float)new_frontier->count/ graph->num_nodes > 0.00001)){
                flag = 0;
                vertex_set_top_to_bottom(new_frontier);
            }
            total_iterations++;
            top_down_iterations++;
        }
        else {
            vertex_set_clear_bottom(new_frontier); 
            bottom_up_step_2(graph, frontier, new_frontier, sol->distances);
            total_iterations++;
        }
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;

    }
    printf("total iterations: %d\n", total_iterations);
    printf("top down iterations: %d\n", top_down_iterations);
    }

// void bfs_hybrid(Graph graph, solution* sol)
// {
//     // //top down 
//     bfs_top_down(graph, sol);

//     // // bottom up 
//     // bfs_bottom_up(graph, sol);
    
//     }


