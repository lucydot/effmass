#! /usr/bin/env python3

import questionary
from effmass import inputs, analysis, extrema, outputs
from random import randint

import warnings
warnings.filterwarnings("ignore")

def cli():

	print("Welcome to effmass 2.0.0 \U0001F388")

	ignore, seedname, fermi_level = None, None, None
	random_int = randint(10000,99999)

	DFT_code = questionary.select(
		"Which DFT code have you used to generate the bandstructure?",
		choices=[
			'Vasp',
			'FHI-aims',
			'Castep'
		]).ask()

	pathname = questionary.path(
		"What's the path to your {} output files?".format(DFT_code),
		default = "./",
		only_directories=True
		).ask()

	if DFT_code == 'Vasp':
		ignore = questionary.text(
		"How many k-points should I ignore at the start of the file? (useful for hybrid calculations)",
		default="0"
		).ask()

	if DFT_code == 'Castep':

		seedname = questionary.text(
		"What's the seedname of your Castep files? (e.g. Si.bands and Si.castep have the seedname Si)",
		).ask()
		fermi_level = questionary.text(
		"I will infer the position of the CBM and VBM from the calculated Fermi level." +
		" If you know this value to be incorrect, please input a more accurate value:",
		).ask()

	extrema_search_depth = questionary.text(
		"How far (in eV) from the CBM (VBM) would you like me to search for minima (maxima)?",
		default = "0.05"
		).ask()

	energy_range = questionary.text(
		"What would you like the energy range (in eV) of each segment to be?",
		default = "0.5"
		).ask()

	which_values = questionary.checkbox(
		"Which values would you like me to calculate?",
		choices=[
			questionary.Choice("parabolic m* (least squares)", checked=True),
			questionary.Choice("parabolic m* (finite difference)", checked=True)
		]).ask()   # need to select oe

	save_plot = questionary.confirm(
		"Would you like me to save a plot of the band segments?",
		default=True,
		auto_enter=False
		).ask()

	save_summary = questionary.confirm(
		"Would you like me to save a summary file?",
		default=True,
		auto_enter=False
		).ask()

	settings = inputs.Settings(extrema_search_depth=float(extrema_search_depth), energy_range=float(energy_range))
	print("Reading in data...")

	if DFT_code == "Vasp":       
		data = inputs.DataVasp(pathname+"/OUTCAR",pathname+"/PROCAR", ignore=int(ignore))

	elif DFT_code == "FHI-aims":
		data = inputs.DataAims(pathname)

	else:
		data = inputs.DataCastep(pathname+"/", seedname)
		if fermi_level:
			data.fermi_level = fermi_level
			data.find_cbm_vbm()
    
	print("Finding extrema...")  
	print("Generating segments...")
	segments = extrema.generate_segments(settings, data)

	print("Calculating effective masses...")
	table = outputs.make_table(segments, which_values)
	outputs.print_terminal_table(table)

	if save_plot:
		print("Plotting segments...")
		outputs.plot_segments(data,settings,segments,savefig=True,random_int=random_int)
		print("Plot of segments saved to effmass_{}.png".format(random_int))

	if save_summary: 
		print("Writing summary file...")

		outputs.print_summary_file(
			random_int,
			DFT_code,
			pathname,
			ignore,
			seedname,
			fermi_level,
			extrema_search_depth,
			energy_range,
			table)

		print("Summary file saved as effmass_{}.txt".format(random_int))


