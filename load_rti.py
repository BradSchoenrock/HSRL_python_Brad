import sys
from hsrl.data_stream.rti import Rti
if __name__ == '__main__':
	print """
sample run:

r=Rti('gvhsrl','04-Feb-12 14:00', 2, 1.5, 15, 2.1, 'flight_plots.json')
	"""
	#r=Rti('gvhsrl','14-Feb-12 23:33', 1, 1.5, 15, 2.1, 'flight_plots.json')
	r=Rti('gvhsrl','19-Feb-12 22:00', 0.5, 0.0, 12, display='flight_plots.json')

#particle size
 g=Rti(('mf2hsrl','magkazrge','magpars2S1','magpars2S2','marinemet','spheroid_particle','multiple_scattering'),'17-jul-13 11:05', '17-jul-13 11:10', 0.0, 1.0, mol_norm_alt=0.2,display='precip_plots.json',particle_parameters='spheroid_particle_parameters_default.json',multiple_scattering='multiple_scattering_paramters_default.json')
	
