#! /usr/bin/env python3

"""
A module for plotting and summarising segments information, 
density-of-states information and effective mass analysis.

"""
import matplotlib
matplotlib.use('Agg',warn=False)

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text

from effmass import ev_to_hartree
from effmass.dos import _check_integrated_dos_loaded
from effmass.dos import _check_dos_loaded
from effmass.analysis import _check_poly_order


def plot_segments(Data, Settings, segments):
    """Plots bandstructure overlaid with the DFT-calculated points for each Segment
    instance. Each Segment is labelled with it's direction in reciprocal space
    and index number from the segments argument.

    Args:
        Data (Data): instance of the :class:`Data` class.
        Settings (Settings): instance of the :class:`Settings` class.
        segments (list(Segment)): A list of instances of the :class:`Segment` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The x-axis of the plot is not to scale.
    """
   
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    [ax.plot(range(len(Data.energies[i])), Data.energies[i] - Data.VBM) for i in range(len(Data.energies))]
    points = [ax.scatter(segments[i].kpoint_indices, segments[i].energies - Data.VBM) for i in range(len(segments))]
    texts = [ax.text(segments[i].kpoint_indices[-1], segments[i].energies[-1] - Data.VBM, str(i)+", "+str(np.round(segments[i].direction,3))) for i in range(len(segments))]    
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), autoalign='x', force_text=(0.01,0.025))

    ax.set_ylim([
        -(Settings.extrema_search_depth + Settings.energy_range + 1),
        (Data.CBM - Data.VBM) +
        (Settings.extrema_search_depth + Settings.energy_range + 1)
    ])

    return fig, ax
    

def plot_integrated_dos(Data):
    """Plots integrated density of states (states/unit-cell) against energy
    (eV).

    Args:
        Data (Data): instance of the :class:`Data` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The valence band maximum is set to 0 eV.
    """
    _check_integrated_dos_loaded(Data)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    energy = [x[0] - Data.VBM for x in Data.integrated_dos]
    dos_data = [x[1] for x in Data.integrated_dos]
    ax.plot(energy, dos_data)
    ax.set_xlabel("Energy, eV")
    ax.set_ylabel("Integrated DOS, states / unit cell")
    ax.axvline(0, linestyle="--")
    ax.axvline(Data.CBM - Data.VBM, linestyle="--")

    return fig, ax


def plot_dos(Data):
    """Plots density of states (states/unit-cell) against energy (eV).

    Args:
        Data (Data): instance of the :class:`Data` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The valence band maximum is set to 0 eV.
    """
    _check_dos_loaded(Data)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    energy = [x[0] - Data.VBM for x in Data.dos]
    dos = [x[1] for x in Data.dos]
    ax.plot(energy, dos)
    ax.set_xlabel("Energy, eV")
    ax.set_ylabel("DOS")
    ax.axvline(0, linestyle="--")
    ax.axvline(Data.CBM - Data.VBM, linestyle="--")

    return fig, ax


def print_results(segment, data, settings, polyfit_order=None):
    
    polyfit_order = settings.degree_bandfit if polyfit_order is None else polyfit_order
    _check_poly_order(polyfit_order)


    print(segment.band_type, segment.direction)
    print("3-point finite difference mass is {:.2f}".format(
        segment.finite_difference_effmass()))
    print("3-point parabolic mass is {:.2f}".format(
        segment.five_point_leastsq_effmass()))
    try: 
        print("weighted parabolic mass is {:.2f}".format(
            segment.weighted_leastsq_effmass()))
    except AssertionError as e:
        print ("----------\n")
        print (e)
        print ("\n-----------")

    try:
        print("alpha is {:.2f} 1/eV".format(
            segment.alpha(polyfit_order=polyfit_order) * ev_to_hartree))
        print("kane mass at bandedge is {:.2f}".format(
            segment.kane_mass_band_edge(polyfit_order=polyfit_order)))
        if segment.explosion_index(polyfit_order=polyfit_order) == len(
                segment.dE_eV):
            print(
                "the Kane quasi linear approximation is valid for the whole segment"
            )
        else:
            print("the Kane quasi-linear approximation is valid until {:.2f} eV".
                  format(segment.dE_eV[segment.explosion_index(
                      polyfit_order=polyfit_order)]))
    except AssertionError as e:
        print ("----------\n")
        print (e)
        print ("\n-----------")

    try:
        print("optical mass at band edge (assuming the Kane dispersion) is {:.2f}".
              format(segment.optical_effmass_kane_dispersion()))
    except AssertionError:
        pass

    plt.figure(figsize=(8, 8))
    plt.plot(
        np.linspace(segment.dk_angs[0], segment.dk_angs[-1], 100),
        np.divide(
            segment.poly_fit(
                polyfit_order=polyfit_order, polyfit_weighting=False),
            ev_to_hartree),
        marker="x",
        ms=5,
        label="polynomial order {}".format(polyfit_order))
    plt.plot(
        np.linspace(segment.dk_angs[0], segment.dk_angs[-1], 100),
        np.divide(segment.finite_difference_fit(), ev_to_hartree),
        marker="<",
        ms=5,
        label="finite diff parabolic")
    plt.plot(
        np.linspace(segment.dk_angs[0], segment.dk_angs[-1], 100),
        np.divide(segment.five_point_leastsq_fit(), ev_to_hartree),
        marker="p",
        ms=5,
        label="five point parabolic")

    try:
        plt.plot(
            np.linspace(segment.dk_angs[0], segment.dk_angs[-1], 100),
            np.divide(segment.weighted_leastsq_fit(), ev_to_hartree),
            marker=">",
            ms=5,
            label="weighted parabolic")
    except AssertionError:
        pass

    try:
        plt.plot(
            np.linspace(segment.dk_angs[0], segment.dk_angs[-1], 100),
            np.divide(
                segment.kane_fit(polyfit_order=polyfit_order), ev_to_hartree),
            marker="o",
            ms=5,
            label="Kane quasi-linear")
    except AssertionError as e:
        pass

    plt.xlabel(r"k ($ \AA^{-1} $)")
    plt.ylabel("energy (eV)")
    plt.scatter(segment.dk_angs, segment.dE_eV, marker="x", s=200, label="DFT")
    plt.legend()
    plt.show()

    fig, axes = plot_segments(data, settings, [segment])
    plt.show()

    plt.figure(figsize=(8, 8))
    idx = segment.explosion_index(polyfit_order=polyfit_order)
    plt.scatter(
        segment.dE_hartree[1:idx + 1],
        segment.transport_effmass(
            polyfit_order=polyfit_order,
            dk=segment.dk_bohr,
            polyfit_weighting=False)[1:idx + 1])
    plt.plot([0, segment.dE_hartree[idx]],
             np.polyval(
                 np.polyfit(
                     segment.dE_hartree[1:idx + 1],
                     segment.transport_effmass(
                         polyfit_order=polyfit_order,
                         dk=segment.dk_bohr,
                         polyfit_weighting=False)[1:idx + 1], 1),
                 [0, segment.dE_hartree[idx]]))
    plt.ylabel("transport mass")
    plt.xlabel("energy (hartree)")
    plt.xlim([0, segment.dE_hartree[idx]])
    plt.show()
