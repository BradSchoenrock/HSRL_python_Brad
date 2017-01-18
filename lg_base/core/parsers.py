#!/usr/bin/env python
import os
import re
from os.path import join

def parse_dates():
    """This module parses directory folder names into dates.
       
    This module works by walking through the directory supplied by HSRL_CONFIG
    and parsing the dates into an array containing sub arrays of year, month and day.
    The 0th index is years, the 1st index is months and the 2nd index is days.
    The month and day sub-indexs contain further sub-indexes which can be 
    accessed by supplying the proper indice which is corresponded by the year. 
    I.e., if you want to access the month data for 2010 and 2010 is index 1
    of years you will access seperated dates as such seperatedDates[1][1]. 
    The first index chooses the month sub-array and the second chooses the year
    2010 in the month index. The same holds true for days but with another index
    for months.
       
    [[Years],[[Months_year0],[Months_year1],[Months_year2]],
    [[[Days_year0,month0],Days_year0,month1...]]]   
    """
    dates = []
    seperatedDates = []
    for root, dirs, files in os.walk(get_workingdirs()):
        for name in dirs:
            f = str(os.path.join(root, name))
            if (f.endswith('raw')):
                dates.append(re.search(
                    '(\d{1,4})[/-](\d{1,2})[/-](\d{1,2})', f).group(0))
    dates.sort()
    split_day = dates.pop(0).split('/')
    seperatedDates.append([split_day[0]])
    seperatedDates.append([[split_day[1]]])
    seperatedDates.append([[[split_day[2]]]])
    i=0
    j=0
    while(len(dates)):
        split_day = dates.pop(0).split('/')
        if (split_day[0] != seperatedDates[0][i]):
            i+=1
            j=0
            seperatedDates[0].append(split_day[0])
            seperatedDates[1].append([split_day[1]])
            seperatedDates[2].append([[split_day[2]]])
        elif (split_day[1] != seperatedDates[1][i][j]):
            j+=1
            seperatedDates[1][i].append(split_day[1])
            seperatedDates[2][i].append([split_day[2]])
        else:
            seperatedDates[2][i][j].append(split_day[2])  
    
    return seperatedDates


def parse_json(json):
    """This module parses the json config files to be imported into enaml.
    
    This module has been made obsolete due to the new format of the json files.
    """
    names=[]
    data=[]
    comments=[]
    labels = []
    commented = False
    labeled = False
    try:
        #file = os.path.join(os.getenv('HSRL_CONFIG'), 'configuration.json')
        f = open(json, 'r')
    except IOError:
        print "configuration.json could not be found in " + \
              os.getenv('HSRL_CONFIG') + "/"
        return [],[],[]
    for line in f:
        line = line.lstrip()
        if re.match('{?"#',line):
            if line.count('#') == 2:
                comment = re.findall('[\w\[\],/:()]+',line)
                comment = " ".join(comment).rstrip(" ,")
                comments.append(comment)
                commented = True
            if line.count('#') == 1:
                label = re.findall('[\w\[\],/:()]+',line)
                label = " ".join(label).rstrip(" ,")
                temp = label.split(':',1)[1].strip('[] ').replace(' ','').split(',')
                for i in range(len(temp)):
                    if temp[i].count('/'):
                        temp[i] = temp[i].split('/')
                labels.append(temp)
                labeled = True                
        if (re.match('{?"(\w+)"', line) != None):
            name = re.search('"(\w+)"', line).group().strip('\"')
            if (re.search('\[(.*?)\]', line) != None):
                value = re.search('\[(.*?)\]', line).group()
                names.append(name)
                data.append(value)
                if not commented:
                    comments.append("")
                if not labeled:
                    labels.append("none")
                commented = False
                labeled = False
    f.close()
    return names, data, comments, labels

def get_workingdirs():
    """This module gets the data directory for given instrument using HRSL_CONFIG.
    
    This module works by opening the directory supplied in the enviroment
    variable HSRL_CONFIG and then parsing it to find the location of the data
    which is needed for the supplied instrument e.g., gvhsrl.
    """
    dir = os.getenv('HSRL_CONFIG')
    instrument = os.getenv('HSRL_INSTRUMENT')
    f = os.path.join(dir, 'hsrl_python.json')
    f = open(f, 'r')
    for line in f:
        if (re.search(instrument, line) != None):
            dir = re.search('(/[\d\w]*)+', line).group()
    return dir    
