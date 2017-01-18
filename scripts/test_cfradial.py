#!/usr/bin/env python
import datetime 
from hsrl.utils.write_cfradial import write_cfradial

start_dt =datetime.datetime(2012, 2, 22, 16, 0, 0)
write_cfradial('cf_radialtest.nc', start_dt, 10)

