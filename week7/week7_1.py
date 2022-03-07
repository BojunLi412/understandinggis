# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 09:45:24 2021

@author: 81465
"""

from networkx import DiGraph, is_eulerian, eulerian_circuit


# Create a graph object
graph = DiGraph()

# Add the 4 islands (nodes), colours as per the above diagram
graph.add_node('Blue')
graph.add_node('Green')
graph.add_node('Yellow')
graph.add_node('Orange')

# Add the 7 bridges (edges), colours as per the above diagram
graph.add_edge('Blue','Green',id=1)	# e.g. this is a bridge from the Blue to Green island
graph.add_edge('Green','Blue',id=2)
graph.add_edge('Green','Orange',id=3)
graph.add_edge('Orange','Yellow',id=4)
graph.add_edge('Yellow','Blue',id=5)
graph.add_edge('Blue','Orange',id=6)
graph.add_edge('Blue','Yellow',id=7)

# Print a report about the graph that you have just made
print(f"We have {len(list(graph.edges()))} bridges between {len(list(graph.nodes()))} islands.")

# Returns True or False to describe whether or not the graph is Eulerian
if is_eulerian(graph) == True:
    
    print("The graph IS Eulerian")
    
    # Return the Eulerian route around Kaliningrad
    print(list(eulerian_circuit(graph)))
    
else:
    print("The graph IS NOT Eulerian")
   



    
    

  
  