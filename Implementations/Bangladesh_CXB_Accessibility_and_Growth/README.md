These notebooks chart the process for analyzing accessibility to services and employment centers in Cox's Bazar, Bangladesh, under various scenarios of investment in transportation infrastructure. The accessibility metrics used are simple (accessibility to the nearest service) and potential (weighted and unweighted accessibility to ALL services). The resulting outputs are used in the forthcoming Cox's Bazar Growth Diagnostic to explore accessibility's relationship to economic growth, job opportunities, and human development in Cox's Bazar.

The notebooks can roughly be broken into 3 parts.

Steps 1 - 3 describe a very standard GOSTNets computation. The only notable difference between code in these notebooks and others in the Implementation folder is the number of destinations, the use of multiple scenarios, and the notable drop in quality in code clarity from GOST's experts to myself.

Steps 4 and 6 further analyze the data. Step 4 prepares population weighted aggregates of the OD matrix outputs at various administrative levels, labels them as necessary, and outputs them into ready-to-use aggregate CSVs and shapefiles. Step 6 does the same but using potential accessibility measures of access, specifically gravity models, in concert with secondary data inputs (prepared separately in QGIS).

Step 5 uses the results of step 4 to relate accessibility aggregates to secondary data in summary charts.

Feel free to direct question to Robert Banick of the ESAPV team (rbanick@worldbank.org)
