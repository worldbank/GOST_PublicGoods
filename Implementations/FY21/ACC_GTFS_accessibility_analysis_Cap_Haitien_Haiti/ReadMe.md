## Phase 1: Accessibility

Notebooks: Use the notebooks the begin with 'Step1...' to 'Step5...'

Traditional approaches had failed to capture the mobility patterns of the poorest and most vulnerable residents who never use transportation services. An innovative approach was needed to address those data gaps. The South-African startup WIMT, supplied an application that was used to draw the first-ever General Transit Feed Specification (GTFS) maps of Tap Taps. Using local surveyors, the application identified stops along the main routes with their associated transit schedules. 

After the main road network is merged with a GTFS tap tap network as well as a network representing ferry travel, accessibility isochrones are created. The number of people are estimated witin a certain travel time of hospitals.

## Phase 2: flooding

Notebooks: Use the notebooks the begin with 'floods...'

Accessibility Analysis During Flood events 

Floods can significantly impact transportation networks and prevent amenities related to education, health care, and employment from being accessible. In order to have a resilient and reliable road system, it is essential to study the vulnerability of the system under several flooding scenarios.  

This project provides a comprehensive analysis of flood impacts on city-scale transportation system topology and accessibility. Accessibility of a road network is evaluated analyzing the ability (alternative routes) to reach amenities under various flood return periods as well as the assessment of the difficulty (an increase of shortest distance) of reaching amenities in the network. 

We plan to run accessibility analysis on a modified graph meant to simulate a flood taking place. Based on data on a 20-year and 50-year horizon we can remove parts of our graph that are projected to be inundated and therefore would impede travel. Using the results, we can compare how much accessibility would be affected by a flood. 

### Centrality Analysis 

Betweenness centrality (or "betweeness centrality") is a measure of centrality in a graph based on shortest paths. For every pair of origins and destinations in a connected graph, there exists at least one shortest path between the vertices such that either the number of edges that the path passes through (for unweighted graphs) or the sum of the weights of the edges (for weighted graphs) is minimized. The betweenness centrality for each edge is the number of these shortest paths that pass through the edge. Edges with a higher edge centrality value represent edges that receive a great share of traffic that passes through them. Using the Python NetworkX library, edge betweenness centrality can be calculated using the edge_betweeness_centrality algorithm. Edge betweenness centrality will be calculated for both the flood scenarios and the non-flood scenarios and compared. Edges that experience a large difference in values would indicate greater changes in the flow of traffic in a flood based scenario and would be useful for transportation planning. Alternative roads could be identified to keep the main economic areas connected despite disruption.  
