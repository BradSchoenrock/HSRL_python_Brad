from hsrl.data_stream.rti import Rti
if __name__ == '__main__':
	#r = Rti('gvhsrl','03-Feb-2012 17:00', 0.25 , 0.0, 15, 2.1, 'all_plots.json')
	#r = Rti('gvhsrl','03-Feb-2012 17:00', 0.25 , 0.0, 15, 2.1, 'flight_plots.json')

	#r=Rti('gvhsrl','13-Feb-12 17:15',0.25 ,0,35.0,4.0,'calibration.json')
	# dfg cal
	#r=Rti('gvhsrl','13-Feb-12 17:00','17:59' ,0,35.0,4.0,display='flight_plots.json')
	#r.cal_gen()
	#     r = Rti('gvhsrl','03-Feb-2012 14:00', 0.25 , 0.0, 15, 2.1, 'flight_plots.json')
	#     r.time_travel('03-Feb-2012 22:15')
	#     r.params('03-Feb-2012 28:00', 0.25, 0.000000, 15.000000, 2.100000 )

	#    r=Rti('gvhsrl','23-jan-12 18:06','18:16',0,15.0,4.0,'all_plots.json')
	#    r = Rti('gvhsrl','21-Jan-12 11:00', 2, 1.5, 15, 2.1,'all_plots.json')
	#    r = Rti('gvhsrl','21-Jan-12 19:30', 0.25, 0, 14.0, 2.0,'flight_plots.json')
	#    r = Rti('gvhsrl','21-Jan-12 12:00', 1, 1.5, 15, 2.1,'all_plots.json')
	#r=Rti('gvhsrl','26-Feb-12 17:00', '17:15', 0.0 , 12.0, 2.1, display='flight_plots.json')
	#r=Rti('gvhsrl','26-Feb-12 17:00', '17:15', 0.0 , 12.0, 2.1, display='all_plots.json')
	#r=Rti('gvhsrl','14-Jul-11 21:12',1 , 0.0 , 12.0, 2.1, display='flight_plots.json')
	r=Rti('gvhsrl','27-Jan-12 19:15', '19:45', 0.0, 8.0, mol_norm_alt=3.1, display='flight_plots.json')
