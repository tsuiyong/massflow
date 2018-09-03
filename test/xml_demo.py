#!/usr/bin/python
# -*- coding:utf-8 -*-

import xml.etree.ElementTree as ET

tree = ET.parse("input.xml")
root = tree.getroot()

for node in root:
    if node.tag == "power_density":
        power_density = node.text
    elif node.tag == "volume":
        volume = float(node.text)
    elif node.tag == "eta":
        eta = float(node.text)
    else:
        Uranium_alias = node.text
