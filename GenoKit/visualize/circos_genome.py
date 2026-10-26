#!/home/zhusitao/anaconda3/bin/python
# -*- coding:utf-8 -*-
'''
filename: circos_genome.py
date: 2025/12/9 19:26
author: Sitao Zhu
mail: zhusitao1990@163.com
'''


import pycircos
import collections
import matplotlib.pyplot as plt


Garc = pycircos.Garc
Gcircle = pycircos.Gcircle
# Set chromosomes
circle = Gcircle()


# pip install python-circos
# https://github.com/ponnhide/pyCircos?tab=readme-ov-file



def chrom_sector(circle, length_path=""):
    # chromosome sector
    # pycircos/example_data_chromosome_general.csv
    with open(length_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            length = int(array[-1])
            # raxis_range: sector position
            arc = Garc(arc_id=name, size=length, interspace=3,
                       raxis_range=(950, 1000), labelposition=60,
                       label_visible=True)
            circle.add_garc(arc)
        circle.set_garcs()



def cytoband(circle, cytoband_path=""):
    # cytoband sector
    color_dict = {"gneg": "#FFFFFF00",
                  "gpos25": "#EEEEEE",
                  "gpos50": "#BBBBBB",
                  "gpos75": "#777777",
                  "gpos100": "#000000",
                  "gvar": "#FFFFFF00",
                  "stalk": "#C01E27",
                  "acen": "#D82322"}
    arcdata_dict = collections.defaultdict(dict)
    # pycircos/example_data_chromosome_cytoband.csv
    with open(cytoband_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            start = int(array[1]) - 1
            width = int(array[2]) - (int(array[1]) - 1)
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["widths"] = []
                arcdata_dict[name]["colors"] = []
            arcdata_dict[name]["positions"].append(start)
            arcdata_dict[name]["widths"].append(width)
            arcdata_dict[name]["colors"].append(color_dict[array[-1]])

    # raxis_range: sector position
    for key in arcdata_dict:
        circle.barplot(key, data=[1] * len(arcdata_dict[key]["positions"]),
                       positions=arcdata_dict[key]["positions"],
                       width=arcdata_dict[key]["widths"],
                       raxis_range=[950, 1000],
                       facecolor=arcdata_dict[key]["colors"])

def scatter(circle, scatter_path="", raxis_range=(860, 940)):
    # scatter plot sector
    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    # pycircos/example_data_point.csv
    with open(scatter_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            start = int(array[1]) - 1
            end = int(array[2])
            mid = (start + end) / 2
            value = float(array[-1])
            values_all.append(value)
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["values"] = []
            arcdata_dict[name]["positions"].append(mid)
            arcdata_dict[name]["values"].append(value)

    vmin, vmax = min(values_all), max(values_all)
    # raxis_range: sector position
    for key in arcdata_dict:
        circle.scatterplot(key, data=arcdata_dict[key]["values"],
                           positions=arcdata_dict[key]["positions"],
                           rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)],
                           raxis_range=raxis_range,
                           facecolor="orangered", spine=True)

def line(circle, line_path="", raxis_range=(770, 850)):
    # line plot sector
    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    # pycircos/example_data_point.csv
    with open(line_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            start = int(array[1]) - 1
            end = int(array[2])
            mid = (start + end) / 2
            value = float(array[-1])
            values_all.append(value)
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["values"] = []
            arcdata_dict[name]["positions"].append(mid)
            arcdata_dict[name]["values"].append(value)

    vmin, vmax = min(values_all), max(values_all)
    # raxis_range: sector position
    for key in arcdata_dict:
        circle.lineplot(key, data=arcdata_dict[key]["values"],
                        positions=arcdata_dict[key]["positions"],
                        rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)],
                        raxis_range=raxis_range,
                        linecolor="royalblue", spine=False)

def bar(circle, bar_path="", raxis_range=(680, 760)):
    # bar plot sector
    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    # "pycircos/example_data_barplot.csv"
    with open(bar_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            start = int(array[1]) - 1
            end = int(array[2])
            width = end - start
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["widths"] = []
                arcdata_dict[name]["values"] = []
            arcdata_dict[name]["positions"].append(start)
            arcdata_dict[name]["widths"].append(width)
            arcdata_dict[name]["values"].append(float(array[-1]))
            values_all.append(float(array[-1]))

    vmin, vmax = min(values_all), max(values_all)
    # raxis_range: sector position
    for key in arcdata_dict:
        circle.barplot(key, data=arcdata_dict[key]["values"],
                       positions=arcdata_dict[key]["positions"],
                       width=arcdata_dict[key]["widths"], base_value=0.0,
                       rlim=[vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)],
                       raxis_range=raxis_range, facecolor="y", spine=True)

def heatmap(circle, heatmap_path="", raxis_range=(630, 670)):
    # heatmap sector
    values_all = []
    arcdata_dict = collections.defaultdict(dict)
    # pycircos/example_data_rect_gradual.csv
    with open(heatmap_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name = array[0]
            start = int(array[1]) - 1
            end = int(array[2])
            width = end - start
            if name not in arcdata_dict:
                arcdata_dict[name]["positions"] = []
                arcdata_dict[name]["widths"] = []
                arcdata_dict[name]["values"] = []
            arcdata_dict[name]["positions"].append(start)
            arcdata_dict[name]["widths"].append(width)
            arcdata_dict[name]["values"].append(float(array[-1]))
            values_all.append(float(array[-1]))

    vmin, vmax = min(values_all), max(values_all)
    # raxis_range: sector position
    for key in arcdata_dict:
        circle.heatmap(key, data=arcdata_dict[key]["values"],
                       positions=arcdata_dict[key]["positions"],
                       width=arcdata_dict[key]["widths"],
                       raxis_range=raxis_range, vmin=vmin, vmax=vmax,
                       cmap=plt.cm.viridis)

def link(circle, link_path="", radius=630):
    # link plot sector
    # pycircos/example_data_links.csv
    with open(link_path, 'r') as f:
        f.readline()
        for line in f:
            array = line.rstrip().split(",")
            name1 = array[0]
            start1 = int(array[1]) - 1
            end1 = int(array[2])
            name2 = array[3]
            start2 = int(array[4]) - 1
            end2 = int(array[5])
            # sector position
            source = (name1, start1, end1, radius)
            destination = (name2, start2, end2, radius)
            circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)
