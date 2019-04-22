# Optimization at World Bank

## Using PuLP to solve Location Science problems


## Initial Location Science Problems to Solve


### P-Median Problem

Objective: Place P facilities to minimize the average travel distance (or time) among all demand points and facilities

-Optionally it can take into account the level of demand at each point (e.g. number of people)

source:

- (Hakimi 1964,  1965)



### Max Cover

Objective: Determine the location of P facilities in order to maximize the demand covered within a pre-specified maximum distance coverage

-Optionally it can take into account the level of demand at each point (e.g. number of people)

source:

- THE MAXIMAL COVERING LOCATION PROBLEM (Church and ReVelle) (in docs)



### Location Set-Covering Problem

Objective: Determine the minimum number of facilities and their locations in order to cover all demands within a pre-specified maximum distance (or time) coverage

source:

- Location Modeling Utilizing Maximum Service Distance Criteria (Church and Meadows) (in docs)


### Partial Set-Covering Problem

Objective: Determine the minimum number of facilities and their locations in order to cover a given fraction of the population within a pre-specified maximum distance (or time) coverage

-Optionally it can take into account the level of demand at each point (e.g. number of people)

source:

- Two New Location Covering Problems: The Partial P-Center Problem and the Partial Set Covering Problem (Daskin and Owen) (in docs)




## Potential Additional Location Science Problems to Solve


### P-Center Problem

Objective: Minimize the maximum travel distance (or time) among all demand points and facilities

-Optionally it can take into account the level of demand at each point (e.g. number of people)

source:

- (Hakimi 1964,  1965)


### Capacitated P-Median Problem /  Fixed charge location problem (FCLP)

Objective: Minimize total facility and transportation costs. In so doing, it determines the optimal number and locations of facilities, as well as the assignments of demand to a facility. Given the fact that the facilities have capacities, demand may not be assigned to its closest facility.

-Optionally it can take into account the level of demand at each point (e.g. number of people)

source:

- Look in Ch 3 of Discrete Network Location Models (in docs) for the Balinski (1965) formulation



### Capacitated Maximal Covering Location Problem (CMCLP) 

Objective: Taking into account the capacity of each facility, determine the location of P facilities in order to maximize the demand covered within a pre-specified maximum distance coverage

-Optionally it can take into account the level of demand at each point (e.g. number of people)


### Location Set-Covering Problem over Time

- more research needed...

source:

- Applications of the Location Set-covering Problem (Revelle and Falkson) (in docs)



### Good PuLP resources:

- http://benalexkeen.com/linear-programming-with-python-and-pulp-part-5/

- An Introduction to pulp for Python Programmers (in docs)



## Limitations of Problem Size due to Combinatorial Complexity

The most important factor affecting how computing time increases as as the size of this problem increases is the number of candidate facilities.

The problems in this repo are being solved using Integer Programming, so there is a limit on the size of the problems that can be solved. This limit of candidate facilities is likely in the hundreds.

To solve larger problems, heuristics should be looked at. For more information, read 'A Computational Comparison of Different Algorithms for Very Large p-median Problems' (in docs)








