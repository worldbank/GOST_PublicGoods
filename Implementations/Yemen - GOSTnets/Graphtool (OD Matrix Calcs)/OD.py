
# coding: utf-8
import sys,os
import networkx as nx
import math
import numpy as np
import time
import pandas as pd
import multiprocessing

import graph_tool as gt
from graph_tool.topology import *
import GOSTnet_graphtool_utils as gn_utils
gt_available = True
print('GT imported')
pth = r'/home/user/home'

sys.path.append(pth)
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def Calculate_OD(G, origins, destinations, fail_value, thread_no, engine, parameter):

    print('\t\t\nThread %d beginning origin-destination calculations' % thread_no)

    OD = np.zeros((len(origins), len(destinations)))
    i = 1

    if engine == 'gt':
        method = 'A2ALL'

        if method == 'A2ALL':
            print('length of thread %s origins set: %s' % (thread_no, len(origins)))
            print('length of thread %s destinations set: %s' % (thread_no, len(destinations)))
            for o in range(0, len(origins)):
                origin = origins[o]
                with time_limit(2):
                    if parameter == 'time_adj':
                        shortest_time = shortest_distance(G, source=origin, weights=G.edge_properties.time_adj)
                    elif parameter == 'walk_time':
                        shortest_time = shortest_distance(G, source=origin, weights=G.edge_properties.walk_time)
                    else:
                        raise ValueError('unrecognized edge parameter! Aborting!')
                        break

                if (o % 50) == 0:
                    print('Reporting: thread %d, calc number %d of %d (%d percent done)' % (
                    thread_no,
                    o,
                    len(origins),
                    (o / len(origins))*100))

                arr = shortest_time.get_array()
                OD[o] = arr[destinations]
                i+=1

    OD_df = pd.DataFrame(OD, columns = destinations, index = origins)

    print('\t\t\nCompleted all calculations for thread %d. length of DF return: %d' % (thread_no, len(OD_df)))

    return OD_df

def Main(passed_dict):
    ### runs all sub functions - generates altered network, simulates journeys,
    ### and then summarises results. Returns a results dictionary for each scenario run

    G = passed_dict['graph']
    origins = passed_dict['origins']
    destinations = passed_dict['destinations']
    fail_value = passed_dict['fail_value']
    thread_no = passed_dict['thread_number']
    engine = passed_dict['engine']
    passemeter = passed_dict['parameter']

    OD = Calculate_OD(G, origins, destinations, fail_value, thread_no, engine, parameter)

    return OD

# ### Section 6: Run Main Function

if __name__ == '__main__':

    combos = [{'graph':'final_G.pickle',
               'parameter':'time_adj',
               'savename':'OD_matrix_driving'},

               {'graph':'final_G_flood.pickle',
                'parameter':'time_adj',
                'savename':'OD_matrix_driving_flood'},

              {'graph':'final_G_walking.pickle',
               'parameter':'walk_time',
               'savename':'OD_matrix_walking'},

              {'graph':'final_G_walking_flood.pickle',
               'parameter':'walk_time',
               'savename':'OD_matrix_walking_flood'}]

    for combo in combos:

        parameter = combo['parameter']
        savename = combo['savename']
        gfile = combo['graph']
        ofile = r'origins_100m_snapped.csv'
        dfiles = [r'schools_snapped.csv', 'health_centers_snapped.csv']

        # #### Add city name to nearest nodes

        G = nx.read_gpickle(os.path.join(pth, gfile))
        print('\nLoaded nx file. Number of nodes: %s. Number of edges: %s' % (G.number_of_nodes(), G.number_of_edges()))

        if gt_available == True:
            print('\ngt available')
            G_gt = gn_utils.nx2gt(G)
            print('\nConverted to gt')
            calc_engine = 'gt'
        else:
            print('\ngt unavailable')
            calc_engine = 'nx'

        print('\nusing %s as calc engine' % calc_engine)

        fail_value = None

        origins_df = pd.read_csv(os.path.join(pth, ofile))
        origins_df['NN'] = origins_df['NN'].astype(int)

        print('\nOrigins loaded. Number of origins: %d' % len(origins_df))

        origins = list(set(list(origins_df.NN)))

        print('\nUnique number of origins: %d' % len(origins))

        unique_ds = []

        for dfile in dfiles:
            destinations_df = pd.read_csv(os.path.join(pth, dfile))
            destinations_df['NN'] = destinations_df['NN'].astype(int)
            destinations = list(set(list(destinations_df.NN)))
            unique_ds.append(destinations)

        destinations = list(set([item for sublist in unique_ds for item in sublist]))

        print('\nDestinations loaded. Number of unique destinations: %d' % len(destinations))

        print('\nUsing how many destinations: %d' % len(destinations))

        ### CHECK ALL NODES ARE REAL ###

        all_node_targets = [origins, destinations]
        unique_nodes = list(set([item for sublist in all_node_targets for item in sublist]))

        for n in unique_nodes:
            if n not in list(G.nodes()):
                print('\nWARNING! Node %d does not exist in the graph. Terminating program immediately' % n)
                sys.exit()
        print('\nCheck complete - all origin / destination nodes exist in graph')

        results = []

        #threads = multiprocessing.cpu_count()
        threads = 10

        pool = multiprocessing.Pool(processes=threads) # start worker processes

        start = time.time()

        print('\nCommencing Main with %s thread(s). Start time: %s' % (threads, time.ctime()))

        # FLIP ORIGINS AND DESTINATIONS FOR SPEED
        if len(origins) > len(destinations):
            flip = 1
            o_2 = destinations
            destinations = origins
            origins = o_2

        print('largest origin: %s' % max(origins))
        print('largest destination: %s' % max(destinations))

        def chunks(l, n):
            for i in range(0, len(l), n):
                yield l[i:i + n]

        d = []

        origins_split = list(chunks(origins, (int(len(origins)/threads)+1)))

        for i in range(0,len(origins_split)):
            if calc_engine == 'nx':
                F = G.copy()
            elif calc_engine == 'gt':
                F = G_gt.copy()

            l = {'graph':F,
                  'origins':origins_split[i],
                  'destinations':destinations.copy(),
                  'fail_value':fail_value,
                  'thread_number':i+1,
                  'engine':calc_engine,
                  'parameter':parameter}

            d.append(l)

        results = pool.map(Main, d)

        elapsed = (time.time() - start)

        print("\t\nAll threads complete after: %s" % (elapsed))

        l = pd.concat(results)

        if flip == 1:
            l = l.transpose()

        l.to_csv(os.path.join(pth, r'{}.csv'.format(savename)))
