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
	"Where's the path to the {} output files?".format(DFT_code),
	default = "./").ask()

print(pathname)