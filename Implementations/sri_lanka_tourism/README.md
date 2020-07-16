#######

# Sri Lanka Transport
- Authors: 
	- Thomas Gertin, Consultant GOST (tgertin@worldbank.org)
	- Walker Kosmidou-Bradley, Geographer ESAPV (wkosmidoubradley@worldbank.org)
	- Johanna Belanger, Consultant ESAPV (jbelanger@worldbank.org)
- Date created: 09.07.2020
- Date updated: 09.07.2020

# Summary of purpose 
This project anaylzes accessibility of tourism destinations to airports and major tourism cities in Sri Lanka. 
The project uses GOSTnets, an open source python wrapper for flexible and large scale network analysis. 
The project analyzes access and provides recomendations for highest value for money of potential road rehabilitation projects.

# Data Catalog

- step1_creating_graph.ipynb
	- creating the graph network from OSM
- step2_assigning_speed.ipynb
	- assigning speeds to the graph network using mapbox telemetry data means and OSM maxspeed attributes by road class.
- step3_calculate_OD.ipynb
	- calculate origin-destination matrices from Bandarainake International Airport to major tourism cities, then from major cities to tourism sites.
	- simulate road improvements using criticality measurement
	- sorts and extracts optimization results using cost budget dictionary. 

# Software
Analysis: GOSTnets (https://github.com/worldbank/GOSTnets - Branch: mapbox-traffic), Jupyter notebooks

# Credits & Acknowledgements
Mapbox 
World Bank Data collaboratives
OpenStreetMap user community

Yeon Soo Kim (Sri Lanka Poverty Team)
Adja Mansora Dahourou, Arnab Bandyopadhyay, Ashini Samarasinghe, Wei Winnie Wang (Sri Lanka Transport Team)
Charles Fox, Benjamin Stewart, Andres Chamorro Elizondo (GOST)

# Licensing & Use limitations
