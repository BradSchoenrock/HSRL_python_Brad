import dplkit.role.decorator
import dplkit.role.narrator
from datetime import datetime,timedelta

@dplkit.role.decorator.exposes_attrs_in_chain(['constants_first','timeinfo'])
@dplkit.role.decorator.autoprovidenested(nestedclasses=[dict],reuseGenerator=False)
class dpl_calvals_narrator(dplkit.role.narrator.aNarrator):
    """ Constants Framestream Narrator

        :param timeinfo: time window info dictionary
        :param calvals: calvals object or a dictionary to use as perpetual constants

        exposed attributes:

        - constants_first (first time's calvals constants)
    """

    @property
    def constants_first(self):
        return self.rs_constants_first

    @property
    def timeinfo(self):
        return self._timeinfo

    def __init__(self,timeinfo,calvals,edgepadding=timedelta(seconds=0)):
        self.edgepadding=edgepadding
        self._timeinfo=timeinfo
        self.calvals=calvals
        if isinstance(calvals,dict):
            self.rs_constants_first=calvals
        else:
            self.rs_constants_first=self.calvals.select_time(self.timeinfo['starttime'])

    def read(self):
        """generator function.
        """
        interval_start_time=self.timeinfo['starttime']-self.edgepadding
        interval_end_time=self.timeinfo['endtime']+self.edgepadding
        
        chunk_start_time=interval_start_time

        now=None
        reverse_padding=self.timeinfo['reverse_padding'] if 'reverse_padding' in self.timeinfo else timedelta(seconds=0)

        #use_end_time= interval_end_time if interval_end_time else (now-reverse_padding)
        while interval_end_time is None or chunk_start_time<interval_end_time:
            now=datetime.utcnow()-reverse_padding
            if isinstance(self.calvals,dict):
                rs_constants = self.calvals.copy()
                rs_constants['next_cal_time']=interval_end_time or now
            else:
                rs_constants = self.calvals.select_time(chunk_start_time)
            if interval_end_time is not None:
                chunk_end_time=min([rs_constants['next_cal_time'],interval_end_time])
            else:
                chunk_end_time=min([rs_constants['next_cal_time'],now])
                if chunk_end_time is now:
                    time.sleep(1)
            if chunk_end_time<=chunk_start_time:
                print 'WARNING dpl_calvals_narrator trying to use 0-length window. going to end of time'
                chunk_end_time=interval_end_time or now
                if chunk_end_time<=chunk_start_time:
                    continue
            yield { 'chunk_start_time':chunk_start_time, 'chunk_end_time':chunk_end_time,
                    'rs_constants':rs_constants}#, 'rs_soundings':rs_soundings}
            #if interval_end_time==None and chunk_end_time==now:
            #    print('bumping into now. Pausing')
            #    sleep(reverse_padding.total_seconds())
            chunk_start_time=chunk_end_time  
        return
