"""
conf.py
Date: 2017/12/28
Author: Jung Nicolas
Description: A simple library to get a configuration object from an ini file
"""

#Native imports
import os

def getconf(conf_file=None):
	if conf_file is None: conf_file = os.path.join(os.path.dirname(__file__),"conf.ini")

	config=dict()
	if not is_ini_formated(conf_file):
		return None
	with open(conf_file) as cfile:
		for line in cfile:
			if line.strip() != "" and line[0] != "#" and line[0] != ";" and line[0] != "[":
				values=line.split("=")
				key=values[0]
				value='='.join(values[1:])
				#config[key.strip()]=value.strip()
				if len(value.strip().split(",")) == 1:
					config[key.strip()]=value.strip()
				elif len(value.strip().split(",")) > 1:
					config[key.strip()]=value.strip().split(",")
	return config

def is_ini_formated(ini_file):
	valid=True
	keywords=[]
	with open(ini_file) as cfile:
		for line in cfile:
			if line.strip() != "" and line[0] != "#" and line[0] != ";" and line[0] != "[":
				parts=line.split("=")
				if len(parts)<2:
					return False
				if parts[0] in keywords:
					return False
				keywords.append(parts[0])
	return True
