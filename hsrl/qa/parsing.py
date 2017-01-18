
#this is where the parsing functions reside.  given a filename, this will convert a file into a dictionary of sections as ordered lists by time (for log entries and bits)
# another function will take the contents of the dictionaries, which are strings describing flags, and convert them into range and/or altitude arrays

from collections import OrderedDict,namedtuple
import numpy as np
import re
from bottleneck import nanmean
from datetime import datetime
import math
import os

BitField=namedtuple('BitField','value mask')

defaultflagbits=dict(Ext=0,BS=1,dep=2) #order of the sets in the array's value. this orders them in both binary and decimal representations

defaultstates=dict(Bad=2,NC=0,Good=1,Caution=4) #bit mask values, as input. these are set according to QA files, and may be merged with other records

reversestates=['NC','Good','Bad','Bad','Caution','Caution','Bad','Bad'] #the states might get merged together, so reversestates must reverse map, setting priorities.

defaultenumerations=dict(NC=0,Good=1,Bad=2,Caution=3) #enumerations of bits, used to make post-merged base-10 values

class fileParser(object):
    def __init__(self):
        self.flagheader=re.compile('^[0-9]{8}, [0-9]{4},')
        self.logheader=re.compile('^[0-9]{8}, [0-9]{4}')
        self.lineno=0
        self.haveLog=False

    def parseDate(self,val):
        numvals=val.split(',')
        try:
            date=int(numvals[0].strip())
            time=int(numvals[1].strip())
        except ValueError:
            raise ValueError('Error parsing "'+val+'" as date and time on line '+ repr(self.lineno))
        return datetime(year=date/10000,month=(date/100)%100,day=date%100,hour=time/100,minute=time%100)

    def file_flagstart(self,basetime):
        return OrderedDict((('header',OrderedDict(date=basetime)),
                            ('content',[])))

    def parseflagheader(self,val):
        return OrderedDict(date=self.parseDate(val))

    def parselogheader(self,val):
        return OrderedDict(date=self.parseDate(val))

    def preparseFlags(self,val):
        print 'Parsing flags "'+val+'"'
        ret=[]
        for _flag in val.strip().split(';'):
            f=OrderedDict()
            flag=_flag.strip()
            if len(flag)==0:
                continue
            print 'Parsing flag "'+flag+'"'
            content=flag.split(',')
            if len(content)!=2 and len(content)!=1:
                raise ValueError('Bad Flag value "'+flag+'" on line '+repr(self.lineno))
            flagtype=content[0].strip().split(' ')
            if len(flagtype)!=2:
                raise ValueError('Bad Flag type "'+content[0]+'" on line '+repr(self.lineno))
            f['name']=flagtype[0]
            f['value']=flagtype[1]
            if len(content)==1:
                f['rangetype']='all'
            else:
                flagrange=content[1].strip().split(' ')
                if len(flagrange)!=2 and len(flagrange)!=3:
                    raise ValueError('Bad Flag range values in "'+flag+'" on line '+repr(self.lineno))
                if len(flagrange)>2:
                    f['rangetype']=flagrange[0]
                    del flagrange[0]
                else:
                    f['rangetype']='msl'#FIXME
                try:
                    f['min']=int(flagrange[0])
                    f['max']=int(flagrange[1])
                except ValueError:
                    raise ValueError('Bad Flag range numbers in "'+flag+'" on line '+repr(self.lineno))
            ret.append(f)
        #FIXME check for full set of flags, full range?
        return ret

    def verifyLastTwo(self,arr):
        if len(arr)<2:
            return
        if arr[-2]['header']['date']>arr[-1]['header']['date']:
            raise RuntimeError('Date on line '+repr(self.lineno)+' is '+repr(arr[-1]['header']['date'])+' which comes before '+repr(arr[-2]['header']['date']))
        if arr[-2]['header']['date']==arr[-1]['header']['date']:
            if len(arr[-2]['content'])==0:
                del arr[-2]
            else:
                RuntimeError('Date on line '+repr(self.lineno)+' is '+repr(arr[-1]['header']['date'])+' which equals prior entry '+repr(arr[-2]['header']['date']))

    def parseLine(self,line,store):
        commentstart=line.find('#')
        if commentstart>=0:
            trimmed=line[:commentstart]
        else:
            trimmed=line
        if len(trimmed.strip())==0:
            return
        if self.flagheader.match(line)!=None:
            m=self.flagheader.match(line)
            store['flags'].append(OrderedDict((('header',self.parseflagheader(line[m.start():m.end()])),('content',self.preparseFlags(line[m.end():].strip())))))
            self.verifyLastTwo(store['flags'])
            self.haveLog=False
        elif self.logheader.match(line)!=None:
            m=self.logheader.match(line)
            store['logs'].append(OrderedDict((('header',self.parselogheader(line[m.start():m.end()])),('content',''))))
            self.verifyLastTwo(store['logs'])
            self.haveLog=True
        else:
            tmp=line.strip()
            if tmp[:3]==';;;':
                return
            if len(tmp)>0:
                if self.haveLog:
                    store['logs'][-1]['content']+=tmp+'\n'
                else:
                    print 'WARNING Unused content on line '+repr(self.lineno)+': "'+tmp+'"'

    def parseFile(self,filename,basetime=None,destination=None):
        ret=destination or OrderedDict((('logs',[]),('flags',[])))
        if basetime:
            ret['flags'].append(self.file_flagstart(basetime))
            self.verifyLastTwo(ret['flags'])
        if not os.access(filename,os.R_OK):
            return ret
        f=file(filename)
        lines=f.readlines()
        for l,val in enumerate(lines):
            self.lineno=l+1
            self.parseLine(val,ret)
        #FIXME do some checks on the resulting dictionary
        #print ret
        return ret

class flagParser(object):
    def __init__(self,altitude_axis,flagbits=None,range_zero_m=None,binwidth=None):
        self.altitude_axis=altitude_axis
        self.altitude_bins=np.zeros([altitude_axis.size+1],dtype=altitude_axis.dtype)
        self.altitude_bins[:-1]=altitude_axis[:]
        self.altitude_bins[-1]=self.altitude_bins[-2]+nanmean(altitude_axis[:-1]-altitude_axis[1:])
        self.binwidth=binwidth#no assumptions
        assert(self.binwidth!=None)
        self.range_zero_m=range_zero_m or 0
        self.flagbits=flagbits or defaultflagbits
        self.flagstates=defaultstates
        self.reversestates=reversestates
        self.flagenumerations=defaultenumerations
        self.state_bitwidth=0
        self.state_bitmask=0
        for f,v in self.flagstates.items():
            self.state_bitmask|=v
        self.state_bitwidth=math.log(self.state_bitmask+1,2)
        self.btecache={}

    @property
    def altitudeAxis(self):
        return self.altitude_axis

    def range_axis(self,maxrange,count=None):
        size=count or ((maxrange+self.range_zero_m)/self.binwidth)
        return np.arange(size)*self.binwidth-self.range_zero_m

    def parseFlags(self,flagdictionary):
        return (self.parseAltitude(flagdictionary),self.parseRange(flagdictionary))#altitude and range

    def flagGen(self,flag):
        m=self.state_bitmask*(2**(self.flagbits[flag['name']]*self.state_bitwidth))
        v=self.flagstates[flag['value']]*(2**(self.flagbits[flag['name']]*self.state_bitwidth))
        return BitField(v,m)

    def parseAltitude(self,flagdictionary):
        ret=np.zeros(self.altitude_axis.shape,dtype='uint32')
        tmp=np.zeros_like(ret)
        axis=self.altitude_bins
        for part in flagdictionary:
            tmp[:]=0
            if part['rangetype'] not in ('msl','all'):
                continue
            val=self.flagGen(part)
            if part['rangetype']=='all':
                tmp[:]=val.value
            else:
                minval=part['min']
                maxval=part['max']
                tmp[np.logical_and(axis[:-1]>=(minval-self.binwidth),axis[1:]<(maxval+self.binwidth))]=val.value
            ret|=tmp
        return ret

    def parseRange(self,flagdictionary):
        ret=None
        for part in flagdictionary:
            if part['rangetype'] not in ('bins','range'):
                continue
            val=self.flagGen(part)
            minval=part['min']
            maxval=part['max']
            if part['rangetype']=='bins':
                minval=minval*self.binwidth-self.range_zero_m
                maxval=maxval*self.binwidth-self.range_zero_m
            tmpr=self.range_axis(maxval)
            if ret==None or ret.size<tmpr.size:
                t=ret
                ret=np.zeros(tmpr.size,dtype='uint32')
                if t!=None:
                    ret[:t.size]=t[:]
            tmp=np.zeros(tmpr.size,dtype='uint32')
            tmp[np.logical_and(tmpr>=(minval-self.binwidth),tmpr<(maxval+self.binwidth))]=val.value
            ret[:tmp.size]|=tmp[:]
        return ret

    def mergeVectors(self,altvals,rangevals,altitude,roll_angle=None,pitch_angle=None):#angles should include the zenith, and corrected for telescope direction (GV)
        if rangevals==None:
            return altvals
        range_axis_as_altitude=self.range_axis(count=rangevals.size)
        if roll_angle!=None and pitch_angle!=None:
            range_axis_as_altitude=range_axis_as_altitude*np.cos(roll_angle) * np.cos(pitch_angle)
        range_axis_as_altitude+=altitude
        ret=altvals.copy()
        for x in self.altitude_axis.size:
            rmask=np.logical_and(range_axis_as_altitude>=self.altitude_bins[x],range_axis_as_altitude<self.altitude_bins[x+1])
            if sum(rmask)==0:
                continue
            ret[x]=np.bitwise_or(altvals[x],np.bitwise_or(rangevals[rmask]))
        return ret

    def convertBitsToEnumeration(self,value):
        #this takes the merged bit field values and replaces it with a base-10 enumeration field (one entry per power of 10)
        if not value in self.btecache:
            print type(value),value
            self.btecache[value]=0
            for n,offset in self.flagbits.items():
                tmp=type(value)(value/(2**(self.state_bitwidth*offset)))&self.state_bitmask
                enumv=self.flagenumerations[self.reversestates[tmp]]
                self.btecache[value]+=enumv*(10**offset)
        return self.btecache[value]

    def translateToEnumeration(self,vector):
        ret=np.zeros_like(vector)
        for x in range(vector.size):
            ret[x]=self.convertBitsToEnumeration(vector[x])
        return ret


def main():
    import sys
    f=fileParser()
    fp=flagParser(np.arange(4000)*7.5,binwidth=7.5)
    x= f.parseFile(sys.argv[1])
    print x
    for f in x['flags']:
        print f['header'],'=',fp.translateToEnumeration(fp.parseFlags(f['content'])[0])
    for f in x['logs']:
        print f['header'],'=',f['content']

if __name__ == '__main__':
    main()