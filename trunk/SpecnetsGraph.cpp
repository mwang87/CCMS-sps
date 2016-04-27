/*
 * SpecnetsGraph.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: Mingxun Wang
 */


#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "FilterableGraph.h"
#include "SpecnetsGraph.h"
#include "utils.h"

using namespace std;


namespace specnets
{
    SpecnetsGraph::SpecnetsGraph()
    : FilterableGraph()
    {
    }
    
    SpecnetsGraph::SpecnetsGraph(SpectrumPairSet pairs_set)
    : FilterableGraph()
    {
        for(int i = 0; i < pairs_set.size(); i++){
            SpectrumPair pair = pairs_set[i];
            int spec1 = pair.spec1;
            int spec2 = pair.spec2;
            
            map<unsigned int, unsigned int>::iterator it;
            
            it=specnets_to_graph_nodes_map.find(spec1);
            
            if(it == specnets_to_graph_nodes_map.end()){
                Node * new_node = addNode();
                specnets_to_graph_nodes_map[spec1] = new_node->getIndex();
                graph_to_specnets_nodes_map[new_node->getIndex()] = spec1;
            }
            
            it=specnets_to_graph_nodes_map.find(spec2);
            
            if(it == specnets_to_graph_nodes_map.end()){
                Node * new_node = addNode();
                specnets_to_graph_nodes_map[spec2] = new_node->getIndex();
                graph_to_specnets_nodes_map[new_node->getIndex()] = spec2;
            }
            
            //Add Edge
            Edge * new_edge = addEdge(specnets_to_graph_nodes_map[spec1], specnets_to_graph_nodes_map[spec2]);
            specnets_to_graph_edge_map[i] = new_edge->getIndex();
            graph_to_specnets_edge_map[new_edge->getIndex()] = i;
        }
        graph_pairs_set = pairs_set;
    }
    
    //Filters out edges in connected component by ever increasing score1 to reduce its size to below maximum
    int SpecnetsGraph::filter_graph_component_size(int maximum_component_size){
        
        bool oversized_component = false;
        
        do
        {
            oversized_component = false;
            
            vector<int> component_sizes;
            vector<unsigned int> components_node;
            
            getComponentSizes(component_sizes, components_node);
            
            for(int i = 0; i < component_sizes.size(); i++){
                if(component_sizes[i] > maximum_component_size){
                    cout<<"COMPONENT TOO BIG AT : "<<component_sizes[i]<<endl;
                    
                    oversized_component = true;
                    
                    vector<Node*> node_list;
                    vector<Edge*> edge_list;
                    
                    clear_traversal_datastructures();
                    getConnectedComponent(getNode(components_node[i]), node_list, edge_list);
                    
                    //cout<<(node_list.size())<<endl;
                    //cout<<(edge_list.size())<<endl;
                    
                    float minimum_component_score = lowest_score_edge(edge_list);
                    
                    cout<<"MINIMUM COMPONENT EDGE SCORE "<<minimum_component_score<<endl;
                    
                    float delta = 0.1;
                    
                    float component_minimum_threshold = minimum_component_score + delta;
                    
                    int removed_edge_count = remove_low_scoring_edges(edge_list, component_minimum_threshold);
                    cout<<"REMOVED " <<removed_edge_count<<" EDGES"<<endl;
                }
            }
            
        }while(oversized_component == true);
        
        return 0;
    }
    
    float SpecnetsGraph::lowest_score_edge(vector<Edge*> edge_list){
        float minimum = 1000.f;
        
        for(int i = 0; i < edge_list.size(); i++){
            //cout<<"edge ptr "<<edge_list[i]<<endl;
            //cout<<"graph edge index "<<edge_list[i]->getIndex()<<endl;
            //cout<<"specnets index "<<graph_to_specnets_edge_map[edge_list[i]->getIndex()]<<endl;
            float score1 = graph_pairs_set[graph_to_specnets_edge_map[edge_list[i]->getIndex()]].score1;
            
            //cout<<i<<"\t"<<score1<<endl;
            
            minimum = min(minimum, score1);
        }
        
        return minimum;
    }
    
    int SpecnetsGraph::remove_low_scoring_edges(vector<Edge*> edge_list, float minimum_threshold){
        
        int number_removed_edges = 0;
        
        for(int i = 0; i < edge_list.size(); i++){
            //cout<<"edge ptr "<<edge_list[i]<<endl;
            //cout<<"graph edge index "<<edge_list[i]->getIndex()<<endl;
            //cout<<"specnets index "<<graph_to_specnets_edge_map[edge_list[i]->getIndex()]<<endl;
            float score1 = graph_pairs_set[graph_to_specnets_edge_map[edge_list[i]->getIndex()]].score1;
            
            if (score1 < minimum_threshold) {
                //Mark this index to be deleted from pairs
                specnets_pairs_deleted.push_back(graph_to_specnets_edge_map[edge_list[i]->getIndex()]);
                
                //Delete Edge From Graph
                removeEdge(edge_list[i]->getIndex());
                number_removed_edges++;
            }
        }
        
        return number_removed_edges;
    }
    
    vector<unsigned int> SpecnetsGraph::get_pairs_deleted(){
        return specnets_pairs_deleted;
    }
    
}