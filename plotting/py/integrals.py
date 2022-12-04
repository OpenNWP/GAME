# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/GAME

# This file is for plotting integrals.

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

fig_save_path = "/home/max/code/GAME/figs"
output_dir = "/home/max/code/GAME/output"
write_out_dry_mass_integral = 1
write_out_rhotheta_integral = 1
write_out_energy_integral = 1
run_id = "ullrich_diff"

# end of usual input section

output_dir = output_dir + "/" + run_id

fig_size = 6
if write_out_dry_mass_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size))
	plt.title("Masses")
	plt.ylabel("change relative to initial value / %")
	data = np.genfromtxt(output_dir + "/masses")
	time_rescale = 1/3600
	time_unit = "hr"
	if (max(data[:, 0]) > 864000):
		time_rescale = 1/86400
		time_unit = "days"
	time_vector = time_rescale*(data[:, 0] - data[0, 0])
	plt.xlabel("time since initialization / " + time_unit)
	plt.xlim([min(time_vector), max(time_vector)])
	data = np.genfromtxt(output_dir + "/masses")
	no_of_constituents = len(data[0, :]) - 1
	if no_of_constituents == 1:
		dry_mass_vector = data[:, 1]
		plt.plot(time_vector, 100.0*(dry_mass_vector/dry_mass_vector[0] - 1.0))
		plt.legend(["Dry mass"])
	if no_of_constituents == 6:
		# dry mass
		dry_mass_vector = data[:, 5] - data[:, 6]
		plt.plot(time_vector, 100.0*(dry_mass_vector/dry_mass_vector[0] - 1.0))
		# the total amount of water in the atmosphere at the beginning
		water_masses_init_sum = data[0, 1] + data[0, 2] + data[0, 3] + data[0, 4] + data[0, 6] + 1.0
		# water vapour
		plt.plot(time_vector, 100.0*(data[:, 6]/(data[0, 6] + 1.0) - 1.0))
		# water in all phases
		plt.plot(time_vector, 100.0*((data[:, 1] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 6])/water_masses_init_sum - 1.0))
		plt.legend(["Dry mass", "Water vapour", "Water (all phases)"])
	plt.grid()
	print("relative dry mass change: " + str(100.0*(dry_mass_vector[-1] - dry_mass_vector[0])/dry_mass_vector[0]) + " %")
	fig.savefig(fig_save_path + "/" + run_id + "_masses_integrals.png", dpi = 500)
	plt.close()
	
if write_out_rhotheta_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size))
	ax = plt.axes()
	ax.grid()
	plt.title("Rho x theta")
	plt.ylabel("change relative to initial value / %")
	data = np.genfromtxt(output_dir + "/potential_temperature_density")
	time_rescale = 1/3600
	time_unit = "hr"
	if (max(data[:, 0]) > 864000):
		time_rescale = 1/86400
		time_unit = "days"
	time_vector = time_rescale*(data[:, 0] - data[0, 0])
	plt.xlabel("time since initialization / " + time_unit)
	plt.xlim([min(time_vector), max(time_vector)])
	entropy_vector = data[:, 1]
	plt.plot(time_vector, 100.0*(entropy_vector/entropy_vector[0] - 1.0))
	fig.savefig(fig_save_path + "/" + run_id + "_rhotheta_integral.png", dpi = 500)
	plt.close()

if write_out_energy_integral == 1:
	fig = plt.figure(figsize = (fig_size, fig_size))
	plt.title("Energy")
	plt.ylabel("change relative to initial value of total energy / %")
	data = np.genfromtxt(output_dir + "/energy")
	time_rescale = 1/3600
	time_unit = "hr"
	if (max(data[:, 0]) > 864000):
		time_rescale = 1/86400
		time_unit = "days"
	time_vector = time_rescale*(data[:, 0] - data[0, 0])
	plt.xlabel("time since initialization / " + time_unit)
	plt.xlim([min(time_vector), max(time_vector)])
	kinetic_vector = data[:, 1]
	potential_vector = data[:, 2]
	internal_vector = data[:, 3]
	total_begin = kinetic_vector[0] + potential_vector[0] + internal_vector[0]
	plt.plot(time_vector, 100.0*(kinetic_vector - kinetic_vector[0])/total_begin)
	plt.plot(time_vector, 100.0*(potential_vector - potential_vector[0])/total_begin)
	plt.plot(time_vector, 100.0*(internal_vector - internal_vector[0])/total_begin)
	plt.plot(time_vector, 100.0*(kinetic_vector + potential_vector + internal_vector - total_begin)/total_begin)
	print("relative energy change: " + str(100.0*(kinetic_vector[-1] + potential_vector[-1] + internal_vector[-1] - total_begin)/total_begin) + " %")
	plt.legend(["kinetic", "potential", "internal", "total"])
	plt.grid()
	fig.savefig(fig_save_path + "/" + run_id + "_energy_integrals.png", dpi = 500)
	plt.close()
	
	
	
	
	
	
	
	
	
	
	

