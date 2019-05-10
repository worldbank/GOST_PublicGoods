def optimize_facility_locations(OD, facilities, p, existing_facilities = None):

    ### Function for identifying spatially optimal locations of facilities ###
    # REQUIRED:   OD - an Origin:Destination matrix, origins as rows, destinations
    #             as columns, in pandas DataFrame format.
    #             facilities - the 'destinations' of the OD-Matrix.
    #             MUST be a list of objects included in OD.columns (or subset)
    #             if certain nodes are unsuitable for facility locations
    #             p - the number of facilities to solve for
    # OPTIONAL:   existing_facilities - facilities to always include in the
    #             solution. MUST be in 'facilities' list
    # -------------------------------------------------------------------------#

    from pulp import LpInteger,LpVariable, LpProblem, lpSum, LpMinimize
    import pandas

    if type(OD) != pandas.core.frame.DataFrame:
        raise ValueError('OD must be pandas Dataframe!')

    for f in facilities:
        if f not in OD.columns:
            raise ValueError('Potential facility locations MUST be in OD.columns')

    if p < 1:
        raise ValueError('need to solve for more than one facility!')
    elif p > len(facilities):
        raise ValueError('need to solve for fewer locations than location options!')

    origins = OD.index
    origins = list(map(int, origins))

    X = LpVariable.dicts('X',(facilities),0,1,LpInteger)

    Y = LpVariable.dicts('Y', (origins,facilities),0,1,LpInteger)

    prob = LpProblem('P Median', LpMinimize)

    prob += sum(sum(OD.loc[i,j] * Y[i][j] for j in facilities) for i in origins)

    prob += lpSum([X[j] for j in facilities]) == p

    for i in origins: prob += sum(Y[i][j] for j in facilities) == 1

    for i in origins:
        for j in facilities:
            prob +=  Y[i][j] <= X[j]

    if existing_facilities != None:
        for e in existing_facilities:
            prob += X[e] == 1

    prob.solve()

    ans = []

    for v in prob.variables():
        subV = v.name.split('_')

        if subV[0] == "X" and v.varValue == 1:
            ans.append(int(str(v).split('_')[1]))

    return ans

def optimize_set_coverage(OD, existing_facilities = None):

    ### Function for identifying spatially optimal locations of facilities ###
    # REQUIRED:   OD - an Origin:Destination matrix, origins as rows, destinations
    #             as columns, in pandas DataFrame format.
    #             facilities - the 'destinations' of the OD-Matrix.
    #             MUST be a list of objects included in OD.columns (or subset)
    #             if certain nodes are unsuitable for facility locations
    # OPTIONAL:   existing_facilities - facilities to always include in the
    #             solution. MUST be in 'facilities' list
    # -------------------------------------------------------------------------#

    from pulp import LpInteger,LpVariable, LpProblem, lpSum, LpMinimize
    import pandas

    origins = OD.index
    origins = list(map(int, origins))

    facilities = OD.keys()
    facilities = list(map(int, facilities))

    X = LpVariable.dicts('X',(facilities),0,1,LpInteger)

    Y = LpVariable.dicts('Y', (origins,facilities),0,1,LpInteger)


    #create a binary variable to state that a facility is placed
    #s = LpVariable.dicts('facility', facilities, lowBound=0,upBound=1,cat=LpInteger)


    prob = LpProblem('Set Cover', LpMinimize)

    #prob += sum(sum(OD.loc[i,j] * Y[i][j] for j in facilities) for i in origins)
    prob += sum(X[j] for j in facilities)


    #for i in origins: prob += sum(Y[i][j] for j in facilities) >= 1

    #find a way to calculate percent coverage

    for i in origins:
        #set of facilities that are eligible to provide coverage to point i
        eligibleFacilities = []
        for j in facilities:
            if OD.loc[i,j] <= 240:
                eligibleFacilities.append(j)
        prob += sum(X[j] for j in eligibleFacilities) >= 1

    prob.solve()

    ans = []

    for v in prob.variables():
        subV = v.name.split('_')

        if subV[0] == "X" and v.varValue == 1:
            ans.append(int(str(v).split('_')[1]))

    return ans
