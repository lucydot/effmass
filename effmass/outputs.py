#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from effmass import dos
from effmass import ev_to_hartree
from effmass.dos import _check_integrated_dos_loaded
from effmass.dos import _check_dos_loaded

plt.style.use('ggplot')
np.set_printoptions(formatter={'float_kind': lambda x: "%.2g" % x})  # print options for all float elements in numpy arrays

def _pretty_info(Segment):
    """
    Helper function which formats segment direction and band data.

    Args:
        Segment (Segment): instance of the :class:`Segment` class.

    Returns:
        str: string containing formated direction and band information.
    """
    info_string =  (
                 "{0:.0f}".format(Segment.direction[0])+
                 "{0:.0f}".format(Segment.direction[1])+
                 "{0:.0f}".format(Segment.direction[2])+"_"+
                 str(Segment.band))
    return info_string

def plot_segments(Data,Settings,segments):
    """
    Plots bandstructure with the DFT-calculated points for each Segment instance overlaid.

    Args:
        Data (Data): instance of the :class:`Data` class.
        Settings (Settings): instance of the :class:`Settings` class.
        segments (list(Segment)): A list of instances of the :class:`Segment` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The x-axis of the plot is not to scale.
    """

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    for i in range(len(Data.energies)):
        ax.plot(range(len(Data.energies[i])),Data.energies[i]-Data.VBM)
    for i in range(len(segments)):
        ax.scatter(segments[i].kpoint_indices,segments[i].energies-Data.VBM)
        ax.annotate(segments[i].direction,xy=(segments[i].kpoint_indices[-1],segments[i].energies[-1]-Data.VBM),xytext=(segments[i].kpoint_indices[-1]+1,segments[i].energies[-1]-Data.VBM+0.1))

    ax.set_ylim([-(Settings.extrema_search_depth+Settings.energy_range+1),(Data.CBM-Data.VBM)+(Settings.extrema_search_depth+Settings.energy_range+1)])

    return fig,ax

def plot_integrated_dos(Data):
    """
    Plots integrated density of states (states/unit-cell) against energy (eV).

    Args:
        Data (Data): instance of the :class:`Data` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The valence band maximum is set to 0 eV.
    """
    _check_integrated_dos_loaded(Data)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    energy = [x[0] - Data.VBM for x in Data.integrated_dos]
    dos_data = [x[1] for x in Data.integrated_dos]
    ax.plot(energy,dos_data)
    ax.set_xlabel("Energy, eV")
    ax.set_ylabel("Integrated DOS, states / unit cell")
    ax.axvline(0,linestyle="--")
    ax.axvline(Data.CBM-Data.VBM,linestyle="--")

    return fig, ax

def plot_dos(Data):
    """
    Plots density of states (states/unit-cell) against energy (eV).

    Args:
        Data (Data): instance of the :class:`Data` class.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.

    Notes:
        The valence band maximum is set to 0 eV.
    """
    _check_dos_loaded(Data)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    energy = [x[0] - Data.VBM for x in Data.dos]
    dos = [x[1] for x in Data.dos]
    ax.plot(energy,dos)
    ax.set_xlabel("Energy, eV")
    ax.set_ylabel("DOS")
    ax.axvline(0,linestyle="--")
    ax.axvline(Data.CBM-Data.VBM,linestyle="--")

    return fig, ax

def plot_fits(Segment,poly_degree,print_to_file=False):
    r"""
    Plots the quadratic fit, polynomial fit and Kane fit to the DFT calculated dispersion.

    The energy is in eV, length in reciprocal space is in angstrom$^{-1}$.

    Args:
        Segment (Segment): instance of the :class:`Segment` class.
        poly_degree (int): the degree of the polynomial fit.
        print_to_file (bool, optional): True to print plot to .png file. Defaults to False.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.
    """
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.scatter(Segment.dk_angs,Segment.dE_eV,marker="x",s=200,label="DFT dispersion")
    ax.plot(np.linspace(Segment.dk_angs[0],Segment.dk_angs[-1],100),[i/ev_to_hartree for i in Segment.return_quadfit(polyfit_order=2,polyfit_weighting=True)],marker="o",label="weighted polynomial parabolic",linewidth=0.01)
    ax.plot(np.linspace(Segment.dk_angs[0],Segment.dk_angs[-1],100),[i/ev_to_hartree for i in Segment.return_finite_difference_fit()],marker="<", ms=5, label="finite difference parabolic")
    ax.plot(np.linspace(Segment.dk_angs[0],Segment.dk_angs[-1],100),[i/ev_to_hartree for i in Segment.return_polyfit(polyfit_order=poly_degree,polyfit_weighting=False)],marker="v",label="polynomial order {}".format(poly_degree))
    if Segment.check_kanefit_points():
        ax.plot(np.linspace(Segment.dk_angs[0],Segment.dk_angs[-1],100), [i/ev_to_hartree for i in Segment.return_kanefit(polyfit_order=poly_degree)], marker="^",label="Kane quasi-linear")
    ax.set_xlim(left=0)

    ax.legend(fontsize=12)

    if print_to_file:
        plt.savefig(_pretty_info(Segment)+"_fit",type="png")

    return fig, ax


def plot_effmass(Segment,poly_degree,print_to_file=False):
    r"""
    Plots the conductivity (second derivative) and transport (first derivative) effective mass against distance from extrema.

    The mass is in units of electron rest mass (:math:`m_e`), distance from extrema is in angstrom :math:`^{-1}`.

    Args:
        Segment (Segment): instance of the :class:`Segment` class.
        Settings (Settings): the degree of the polynomial fit used to calculate the effective masses.
        print_to_file (bool, optional): True to print plot to .png file. Defaults to False.

    Returns:
        Figure, Axes: tuple containing instance of the `matplotlib.pyplot.figure <https://matplotlib.org/api/figure_api.html>`_ class and `matplotlib.pyplot.axes <https://matplotlib.org/api/axes_api.html>`_ class.
    """
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.plot(Segment.dk_angs,Segment.calc_conductivity_effmass(polyfit_order=poly_degree,dk=Segment.dk_angs),label="conductivity mass")
    ax.plot(Segment.dk_angs,Segment.calc_transport_effmass(polyfit_order=poly_degree,dk=Segment.dk_angs),label="transport mass")

    ax.set_xlim(left=0)

    ax.legend(fontsize=12)

    if print_to_file:
        plt.savefig(_pretty_info(Segment)+"_effmass",type="png")    

    return fig, ax

def print_results(segment,data,settings,polyfit_order=6):

    print (segment.ptype, segment.direction)
    print ("finite difference mass is {:.2f}".format(segment.calc_finite_difference_effmass()))
    print ("3-point parabolic mass is {:.2f}".format(segment.calc_five_point_effmass())) 
    print ("weighted parabolic mass is {:.2f}".format(segment.calc_conductivity_effmass(polyfit_order=2,polyfit_weighting=True)[0]))
    if segment.calc_alpha(polyfit_order=polyfit_order) is not None:
        print ("alpha is {:.3f} 1/eV".format(segment.calc_alpha(polyfit_order=polyfit_order)*ev_to_hartree))
        print ("kane mass at bandedge is {:.3f}".format(segment.calc_kane_mass_bandedge(polyfit_order=polyfit_order)))
    else: 
        print ("no kane parameters calculated")
    if segment.explosion_index(polyfit_order=polyfit_order) == len(segment.dE_eV):
        print ("quasi linear approximation valid for whole segment")
        
    else:
        print ("range of quasi-linear until {:.3f} eV".format(segment.dE_eV[segment.explosion_index(polyfit_order=polyfit_order)]))
    if segment.ptype == "electron":
        fermi_level = data.CBM
    if segment.ptype == "hole":
        fermi_level = data.VBM
    print ("optical mass at band edge (assuming Kane dispersion) is {:.2f}".format(segment.calc_optical_effmass_kane_dispersion(fermi_level = fermi_level)))    
    
    plt.figure(figsize=(8,8))
    plt.plot(np.linspace(segment.dk_angs[0],segment.dk_angs[-1],100),np.divide(segment.return_polyfit(polyfit_order=polyfit_order,polyfit_weighting=False),ev_to_hartree),marker="x",ms = 5,label="polynomial order {}".format(polyfit_order))
    plt.plot(np.linspace(segment.dk_angs[0],segment.dk_angs[-1],100),np.divide(segment.return_finite_difference_fit(),ev_to_hartree),marker="<", ms=5, label="finite diff parabolic")
    plt.plot(np.linspace(segment.dk_angs[0],segment.dk_angs[-1],100),np.divide(segment.return_five_point_fit(),ev_to_hartree),marker="p",ms = 5,label="five point parabolic")
    plt.plot(np.linspace(segment.dk_angs[0],segment.dk_angs[-1],100),np.divide(segment.return_polyfit(polyfit_order=2,polyfit_weighting=True),ev_to_hartree),marker=">",ms = 5,label="weighted parabolic")
    if segment.return_kanefit(polyfit_order=polyfit_order) is not None:
        plt.plot(np.linspace(segment.dk_angs[0],segment.dk_angs[-1],100),np.divide(segment.return_kanefit(polyfit_order=polyfit_order),ev_to_hartree),marker="o",ms = 5,label="Kane quasi-linear")
    plt.xlabel(r"k ($ \AA^{-1} $)")
    plt.ylabel("energy (eV)")
    plt.scatter(segment.dk_angs,segment.dE_eV,marker="x", s=200,label="DFT")
    plt.legend()
    plt.show()
    
    fig,axes = plot_segments(data,settings,[segment])
    plt.show()

    plt.figure(figsize=(8,8))
    idx = segment.explosion_index(polyfit_order=polyfit_order)
    plt.scatter(segment.dE_hartree[1:idx+1],segment.calc_transport_effmass(polyfit_order=polyfit_order,dk=segment.dk_bohr,polyfit_weighting=False)[1:idx+1])
    plt.plot([0,segment.dE_hartree[idx]],np.polyval(np.polyfit(segment.dE_hartree[1:idx+1],segment.calc_transport_effmass(polyfit_order=polyfit_order,dk=segment.dk_bohr,polyfit_weighting=False)[1:idx+1],1),[0,segment.dE_hartree[idx]]))
    plt.ylabel("transport mass")
    plt.xlabel("energy (hartree)")
    plt.xlim([0,segment.dE_hartree[idx]])
    plt.show()

    
