
def fmt_mpl_datetime(mpl_date, fmt = '%Y-%m-%d %H:%M:%S'):
    """ format Matplotlib dates into strings """
    from matplotlib.dates import num2date
    tmp=mpl_date
    if type(tmp) is float:
#        print n
        tmp=num2date(mpl_date)
    return tmp.strftime(fmt)

def fmt_mpl_time(mpl_date):
    return fmt_mpl_datetime(mpl_date, fmt = '%H:%M:%S')

