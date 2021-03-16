#! /usr/bin/env python3

import questionary

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

if DFT_code == 'Castep':

	seedname = questionary.text(
	"What's the seedname of your Castep files? (e.g. Si.bands and Si.castep have the seedname Si)",
	).ask()
	fermi_level = questionary.text(
	"I will infer the position of the CBM and VBM from the calculated Fermi level." +
	" If you know this value to be incorrect, please input a more accurate value:",
	).ask()

extrema_search_depth = questionary.text(
	"How far from the CBM (VBM) would you like me to search for minima (maxima)")

energy_range = questionary.text(
	"What would you like the energy range of each segment to be?")

which_values = questionary.checkbox(
	"Which values would you like me to calculate?",
	choices=[
		Choice("parabolic m* (least squares)", checked=True),
		Choice("parabolic m* (finite difference)", checked=True)
	]).ask()

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

search_depth = questionary.text(
	"")

# Output table:   particle band-index direction least-squares m* finite-difference m*

