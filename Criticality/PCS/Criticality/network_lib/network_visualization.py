from __future__ import division
from matplotlib import pyplot as plt
from matplotlib.pylab import *
import matplotlib.colors as colors

import pandas as pd
import math
import numpy as np

import geopandas as gp

__all__ = ['plot_network_admcolmap_betweenness',
           'plot_socioeconomic_attribute',
           'truncate_colormap',
           'plot_network_admcolmap_betweenness_new',
           'plot_od_heatmap']

def plot_network_admcolmap_betweenness(gdf,gdf2, colname,betweenness_string,
                                       cmap='OrRd', linewidth=1.25, edgecolor='grey',
                                       maxbetweenness=0, maxpop=0, thres1=0.1, thres2=0.2):
    fig, ax = plt.subplots(figsize=(12,9))

    ax.set_aspect('equal')

    valmin1 = min(list(gdf2[betweenness_string]))
    valmax1 = max(list(gdf2[betweenness_string]))
    gdf2.plot(ax=ax, column=betweenness_string, cmap=cmap,vmin=valmin1, vmax=valmax1, linewidth=linewidth)

    #adjust linewidth based on betweenness
    betweenness_list = list(gdf2[betweenness_string])
    #change small betweenness values to 0.1 so that they are still visible in the figure
    betweenness_list = [1 if x < thres1 else 2 if x >= thres1 and x < thres2 else 3.5 for x in betweenness_list]
    #betweenness_list = [0.1 if x < 0.1 else x for x in betweenness_list]
    i = 0
    for ln in ax.lines:
        ln.set_linewidth(betweenness_list[i]*1)
        ln.set_linewidth(betweenness_list[i])
        i +=1

    valmin2 = min(list(gdf[colname]))
    valmax2 = max(list(gdf[colname]))
    gdf.plot(ax=ax, column=colname, cmap='Greys',vmin=valmin2, vmax=valmax2, linewidth=0.5, edgecolor=edgecolor, alpha=0.3)

    ax.set_title(colname)
    #remove the lon-lat in the x-y axis of the plot
    ax.axis('off')

    # add colorbar1
    fig = ax.get_figure()
    cax = fig.add_axes([0.85, 0.45, 0.02, 0.43])
    sm = plt.cm.ScalarMappable(cmap='Greys')
    columnlist = list(gdf[colname])
    columnlist.append(0)
    columnlist.append(maxpop) #hardcoded, not good
    cbmin, cbmax = min(columnlist), max(columnlist)
    sm.set_array(columnlist)
    cb = plt.colorbar(sm, cax=cax, label = colname, alpha=0.3)
    labels = [0, cbmax/4, cbmax/4*2, cbmax/4*3, cbmax/4*4]
    loc = labels
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=16))
    cb.ax.tick_params(labelsize=16)

    #add colorbar2
    fig = ax.get_figure()
    cax = fig.add_axes([0.7, 0.45, 0.02, 0.43])
    sm = plt.cm.ScalarMappable(cmap=cmap)
    columnlist = list(gdf2[betweenness_string])
    columnlist.append(0)
    columnlist.append(maxbetweenness)
    cbmin, cbmax = min(columnlist), max(columnlist)
    cbmin, cbmax = round(cbmin,3), round(cbmax,3)
    sm.set_array(columnlist)
    cb = plt.colorbar(sm, cax=cax, label=betweenness_string)
    labels = [0, cbmax/4, cbmax/4*2, cbmax/4*3, cbmax/4*4]
    loc = labels
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=16))
    cb.ax.tick_params(labelsize=16)

def plot_socioeconomic_attribute(gdf, colname,cmap='OrRd', linewidth=1.25, edgecolor='grey', maxpop=0):
    print('maximum number of '+colname+' is',max(list(gdf[colname])))
    fig, ax = plt.subplots(figsize=(12,9))

    ax.set_aspect('equal')

    valmin2 = min(list(gdf[colname]))
    valmax2 = max(list(gdf[colname]))
    gdf.plot(ax=ax, column=colname, cmap=cmap,vmin=valmin2, vmax=valmax2, linewidth=0.5, edgecolor=edgecolor, alpha=0.3)

    ax.set_title(colname)

    #remove the lon-lat in the x-y axis of the plot
    ax.axis('off')

    # add colorbar1
    fig = ax.get_figure()
    cax = fig.add_axes([0.7, 0.45, 0.02, 0.43])
    sm = plt.cm.ScalarMappable(cmap=cmap)
    columnlist = list(gdf[colname])
    columnlist.append(0)
    columnlist.append(maxpop)
    cbmin, cbmax = min(columnlist), max(columnlist)
    sm.set_array(columnlist)
    cb = plt.colorbar(sm, cax=cax, label = colname, alpha=0.3)
    labels = [0, cbmax/4, cbmax/4*2, cbmax/4*3, cbmax/4*4]
    loc = labels
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=16))
    cb.ax.tick_params(labelsize=16)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def _get_percentile(gdf2, col, n):
    #get n-th percentile of a DataFrame column
    get_col = list(gdf2[col])
    get_col = [x for x in get_col if x > 0]
    nth_percentile = np.percentile(get_col, n)

    return nth_percentile

def plot_network_admcolmap_betweenness_new(gdf, gdf2, colname,betweenness_string,
                                       cmap='OrRd', linewidth=1.25, edgecolor='grey',
                                       maxbetweenness=0, maxpop=0, perc1=60, perc2=90):
    fig, ax = plt.subplots(figsize=(12,9))

    ax.set_aspect('equal')

    valmin1 = min(list(gdf2[betweenness_string]))
    valmax1 = max(list(gdf2[betweenness_string]))
    thres1 = _get_percentile(gdf2, betweenness_string, perc1)
    thres2 = _get_percentile(gdf2, betweenness_string, perc2)

    gdf2.plot(ax=ax, column=betweenness_string, cmap=cmap,vmin=valmin1, vmax=valmax1, linewidth=linewidth)

    #adjust linewidth based on betweenness
    betweenness_list = list(gdf2[betweenness_string])
    #change the linewidth based on the percentile
    betweenness_list = [1 if x < thres1 else 2 if x >= thres1 and x < thres2 else 3 for x in betweenness_list]
    i = 0
    for ln in ax.lines:
        ln.set_linewidth(betweenness_list[i]*1)
        i +=1

    valmin2 = min(list(gdf[colname]))
    valmax2 = max(list(gdf[colname]))
    gdf.plot(ax=ax, column=colname, cmap='Greys',vmin=valmin2, vmax=valmax2, linewidth=0.5, edgecolor=edgecolor, alpha=0.3)

    ax.set_title(colname)
    #remove the lon-lat in the x-y axis of the plot
    ax.axis('off')

    # add colorbar1
    fig = ax.get_figure()
    cax = fig.add_axes([0.85, 0.45, 0.02, 0.43])
    sm = plt.cm.ScalarMappable(cmap='Greys')
    columnlist = list(gdf[colname])
    columnlist.append(0)
    columnlist.append(maxpop) #hardcoded, not good
    cbmin, cbmax = min(columnlist), max(columnlist)
    sm.set_array(columnlist)
    cb = plt.colorbar(sm, cax=cax, label = colname, alpha=0.3)
    labels = [0, cbmax/4, cbmax/4*2, cbmax/4*3, cbmax/4*4]
    loc = labels
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=16))
    cb.ax.tick_params(labelsize=16)

    #add colorbar2
    fig = ax.get_figure()
    cax = fig.add_axes([0.7, 0.45, 0.02, 0.43])
    sm = plt.cm.ScalarMappable(cmap=cmap)
    columnlist = list(gdf2[betweenness_string])
#     columnlist.append(0)
    columnlist.append(maxbetweenness)
    cbmin, cbmax = min(columnlist), max(columnlist)
#     cbmin, cbmax = round(cbmin,3), round(cbmax,3)
    sm.set_array(columnlist)
    cb = plt.colorbar(sm, cax=cax, label=betweenness_string)
    poin1 = cbmin+(cbmax-cbmin)/4
    poin2 = cbmin+(cbmax-cbmin)/4*2
    poin3 = cbmin+(cbmax-cbmin)/4*3
    labels = [cbmin, poin1, poin2, poin3, cbmax]
    loc = labels
    cb.set_ticks(loc)
    cb.set_ticklabels(labels)
    cb.ax.yaxis.label.set_font_properties(matplotlib.font_manager.FontProperties(size=16))
    cb.ax.tick_params(labelsize=16)

def _log_base_n(x,logn):
    try:
        return math.log(x,logn)
    except:
        return 0

def plot_od_heatmap(OD_df, gdf_points, log=False, logn=100, division=False):

    #adopted from http://nbviewer.jupyter.org/gist/joelotz/5427209

    #Scale data logarithmically if we want to dampen the Chittagong effect (tremendous amount of goods
    #is transported to Chittagong)
    if log:
        OD_df = OD_df.applymap(lambda x: _log_base_n(x, logn))

    # Plot it out
    fig, ax = plt.subplots()

    #if we don't want to aggregate to division level
    if not division:
        heatmap = ax.pcolor(OD_df, cmap=plt.cm.Blues, alpha=0.8)

        ##################################################
        ## FORMAT ##
        ##################################################

        fig = plt.gcf()
        fig.set_size_inches(14,14)

        # turn off the frame
        ax.set_frame_on(False)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(OD_df.shape[0])+0.5, minor=False)
        ax.set_xticks(np.arange(OD_df.shape[1])+0.5, minor=False)

        # want a more natural, table-like display
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # Set the labels
        ax.set_xticklabels(gdf_points.District, minor=False)
        ax.set_yticklabels(gdf_points.District, minor=False)

    #if we want to aggregate to division level
    else:
        OD_dummy = OD_df.copy()
        gdf_points_dummy = gdf_points.copy()

        node_division_dict = dict(zip(list(gdf_points_dummy['Node']), list(gdf_points_dummy['Division'])))

        OD_dummy.index = [node_division_dict[x] for x in OD_dummy.columns]
        OD_dummy.columns = [node_division_dict[x] for x in OD_dummy.columns]

        OD_dummy = OD_dummy.groupby(OD_dummy.index).sum().groupby(OD_dummy.columns, axis=1).sum()

        heatmap = ax.pcolor(OD_dummy, cmap=plt.cm.Blues, alpha=0.8)

        ##################################################
        ## FORMAT ##
        ##################################################

        fig = plt.gcf()
        fig.set_size_inches(14,14)

        # turn off the frame
        ax.set_frame_on(False)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(OD_dummy.shape[0])+0.5, minor=False)
        ax.set_xticks(np.arange(OD_dummy.shape[1])+0.5, minor=False)

        # want a more natural, table-like display
        ax.invert_yaxis()
        ax.xaxis.tick_top()

        # Set the labels
        ax.set_xticklabels(OD_dummy.columns, minor=False, fontsize=18)
        ax.set_yticklabels(OD_dummy.index, minor=False, fontsize=18)

    # rotate the labels
    plt.xticks(rotation=90)

    # give the x and y label
    plt.xlabel('To', fontsize=18)
    ax.xaxis.set_label_position('top')
    plt.ylabel('From', fontsize=18)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
