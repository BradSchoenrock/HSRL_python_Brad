
dpl_rti.dpl_rti is a wrapper to rti.Rti in the form of a DPL-compatible generator, returning the processed data structure generated, including all intermediate processed quantities

parameters are instrument (e.g. 'bagohsrl'), start and end times of the request (as datetime.datetime objects), time resolution (unused) and max time window returned per iteration (may be less) as datetime.timedelta objects, min and max altitudes in meters, and altitude resolution (unused) in meters
