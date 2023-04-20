#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 14:49:36 2022

@author: egorgusev
"""

from dijkstar import Graph, find_path
import pandas as pd
from tqdm import tqdm
import geopandas as gpd
from shapely.ops import nearest_points
from datetime import timedelta
import movingpandas as mpd
from lxml import etree
from os import path
import warnings
import glob
from pyproj import Transformer
warnings.filterwarnings("ignore", category=FutureWarning)

# %%

utm_crs = 'EPSG:32635'
wgs84 = 'EPSG:4326'
transformer = Transformer.from_crs(utm_crs, wgs84)

def transform_crs(data, transformer):


    lon = transformer.transform(data['x'], data['y'])[0]
    lat = transformer.transform(data['x'], data['y'])[1]

    data['y'] = lon
    data['x'] = lat
    return data

# %%

network = 'network.xml'
network = etree.parse(network)

nodes = {}
for node in network.find('nodes').findall('node'):
    coords = {'x': float(node.get('x')), 'y': float(node.get('y'))}
    nodes[node.get('id')] = coords

nodes = pd.DataFrame(nodes).T
nodes = gpd.GeoDataFrame(nodes, geometry=gpd.points_from_xy(nodes['x'], nodes['y'], crs=utm_crs))

graph = Graph()

links = {}
for link in network.find('links').findall('link'):
    node_from = link.get('from')
    node_to = link.get('to')
    length = float(link.get('length'))
    weight = length/float(link.get('freespeed'))
    links[link.get('id')] = (node_from, node_to, weight)
    graph.add_edge(node_from, node_to, (length, weight))

def cost_func(u, v, edge, prev_edge):
    length, time_sec = edge
    return time_sec

nodes_union = nodes.geometry.unary_union

#%%
date = '2022-08-23'
ods = []
traj_split = []
trajs_length = []
all_files = glob.glob(path.join(date + '/tracks/*.csv'))
pbar = tqdm(len(all_files))
print(len(all_files))

for filename in all_files:
    pbar.update(1)
    try:
        tracks = gpd.read_file(filename, index_col=None, header=0)
    except:
        print('read_file error')
        pass
    tracks['time'] = pd.to_datetime(tracks['time'])
    tracks = tracks[(tracks.time.dt.hour < 6)]

    tracks.geometry = gpd.points_from_xy(tracks['lon'], tracks['lat'], crs = 'EPSG:3857')
    tracks_wgs84 = tracks.to_crs('EPSG:4326')

    tracks_wgs84['lon'] = tracks_wgs84.geometry.x
    tracks_wgs84['lat'] = tracks_wgs84.geometry.y
    
    if len(tracks_wgs84) > 2:
        trajectory = mpd.Trajectory(tracks_wgs84, 'vehicleid', t='time', x='lon', y='lat')
        
        trajectory.add_speed(overwrite=True)
        trajectory.df
        trajectory.df['speed'] = trajectory.df['speed']*3.6

        gap_split = mpd.ObservationGapSplitter(trajectory).split(gap=timedelta(minutes=5))
        for traj in gap_split:      
            split = mpd.StopSplitter(traj).split(max_diameter=300, min_duration=timedelta(minutes=2), min_length=5000)
            if len(split) > 0:
                 origins = split.get_start_locations().reset_index()
                 origins.crs = wgs84
                 origins = origins.to_crs(utm_crs)
                 origins = origins.add_suffix('_origin')
                 dests = split.get_end_locations().reset_index()
                 dests.crs = wgs84
                 dests = dests.to_crs(utm_crs)
                 dests = dests.add_suffix('_dest')
                 clean_orig = []
                 clean_dest = []
                 traj_gdf = split.to_traj_gdf()
                 traj_gdf['vehicleid'] = origins['vehicleid_origin']
                 for index, origin in origins.iterrows():
                     try:
                         dest = dests.loc[index]
                         traj_length = split.to_traj_gdf().loc[index]['length']
                         
                         nearest_from = nodes.geometry == nearest_points(origin['geometry_origin'], nodes_union)[1]
                         nearest_from = nodes[nearest_from].index[0]
                         
                         nearest_to = nodes.geometry == nearest_points(dest['geometry_dest'], nodes_union)[1]
                         nearest_to = nodes[nearest_to].index[0]
                         
                         length_sum = 0
                         for edge in find_path(graph, nearest_from, nearest_to, cost_func=cost_func)[1]:
                             length_sum+=edge[0]
                         if length_sum > 0:
                             if abs(traj_length - length_sum)/length_sum < 0.15:
                                 if (origin['speed_origin'] < 30) and (dest['speed_dest'] < 30):
                                     clean_orig.append(origin)
                                     clean_dest.append(dest)      
                                     traj_split.append(traj_gdf.loc[index])
                     except:
                         print('error')
                         pass
                 od = pd.concat([pd.DataFrame(clean_orig), pd.DataFrame(clean_dest)], axis=1)
                 if len(od) > 0:
                     ods.append(od)
pbar.close()

trajs = pd.concat(traj_split, axis=1).T.drop_duplicates()
ods = pd.concat(ods)
ods['hour_origin'] = pd.to_datetime(ods['index_origin']).dt.hour
ods = ods.drop_duplicates()

trajs.to_csv('trajectories_{}.csv'.format(date))
ods.to_csv('O-D_matrix_{}.csv'.format(date))