# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 08:16:20 2016

@author: GDM
"""

import time
import datetime


def getDay():
    return (datetime.date.today()).isoformat()
    
def getCurrTime():
        currTime = time.strftime('%y-%m-%d %H:%M:%S', time.localtime(time.time()))
        return currTime

def computeRunTime(startTime, endTime):
    start = datetime.datetime.strptime(startTime, '%y-%m-%d %H:%M:%S')
    start_sec_float = time.mktime(start.timetuple())
    
    end = datetime.datetime.strptime(endTime, '%y-%m-%d %H:%M:%S')
    end_sec_float = time.mktime(end.timetuple())

    return end_sec_float - start_sec_float
    
    