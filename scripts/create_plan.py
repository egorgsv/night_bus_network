#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 13:50:00 2022

@author: egorgusev
"""

# %% import
import pandas as pd
from lxml import etree
from pyproj import Transformer
import warnings
warnings.filterwarnings("ignore")
#%%

utm_epsg = 'EPSG:32635'
wgs84_epsg = 'EPSG:4326'
transformer = Transformer.from_crs(
    wgs84_epsg, utm_epsg)  # DO NOT FORGET TO CHANGE CRS

transformer_to_wgs = Transformer.from_crs(
    utm_epsg, wgs84_epsg) 


def wgs84_to_utm(x, y):
    lon, lat = transformer.transform(y, x)
    return (str(lon), str(lat))

def utm_to_wgs84(x, y):
    lat, lon = transformer_to_wgs.transform(x, y)
    return (str(lon), str(lat))


#%%

date = '2022-08-23'
od = pd.read_csv('O-D_matrix_{}.csv'.format(date))

od['time_origin'] = pd.to_datetime(od['index_origin']).dt.strftime('%H:%M')
plans = etree.Element('population')

for i, row in od.iterrows():
    person = etree.SubElement(plans, 'person')
    person.attrib['id'] = str(i)
    plan = etree.SubElement(person, 'plan')
    act1 = etree.SubElement(plan, 'activity')
    act1.attrib['type'] = 'origin'
    lon, lat = wgs84_to_utm(row['lon_origin'], row['lat_origin'])
    act1.attrib['x'] = lon
    act1.attrib['y'] = lat
    act1.attrib['end_time'] = row['time_origin']
    leg = etree.SubElement(plan, 'leg')
    leg.attrib['mode'] = 'pt'
    act2 = etree.SubElement(plan, 'activity')
    lon, lat = wgs84_to_utm(row['lon_dest'], row['lat_dest'])
    act2.attrib['type'] = 'destination'
    act2.attrib['x'] = lon
    act2.attrib['y'] = lat
    

with open('taxi_plans_{}.xml'.format(date), 'wb') as f:
    f.write(etree.tostring(plans, pretty_print=True,
            doctype='<!DOCTYPE population SYSTEM "http://www.matsim.org/files/dtd/population_v6.dtd">'))

