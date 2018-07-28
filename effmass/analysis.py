#! /usr/bin/env python3
"""
A module for analysing the data contained in a :class:`Segment` object.

Contains the :class:`Segment` class and methods for calculating various definitions of the effective mass.
"""

import warnings
import numpy as np
from scipy import integrate
from effmass import extrema
from effmass import angstrom_to_bohr
from effmass import ev_to_hartree
from effmass import boltzmann
from effmass import q

warnings.simplefilter(action='ignore', category=FutureWarning)


def _check_poly_order(polyfit_order):
    """Raises an AssertionError if the order of the polynomial is less than 2.

    Args:
        polyfit_order (int): the order of the polynomial fit.

    Returns:
        None.
    """

    if polyfit_order < 2:
        raise ValueError("the order of the polynomial must be 2 or higher")


def _solve_quadratic(a, b, c):
    r"""
    Solves quadratic equation of the form :math:`ax^2+bx+c=0` for multiple values of c.
    If the determinant is more than 0 (two solutions), it always returns the larger root.

    Args:
        a (float): the coefficient of :math:`x^2`.
        b (float): the coefficient of :math:`x^1`.
        c (list(float)): a list of coefficients of :math:`x^0`.

    Returns:
        list(float): the root of the equation for each value of c.
    
    """
    det = [b**2 - 4 * a * item for item in c]
    x = []
    for d in det:
        if d < 0:
            x.append(np.nan)
        if d == 0:
            x.append(-b / (2 * a))
        if d > 0:
            x.append((-b + np.sqrt(d)) / (2 * a))
    return x


class Segment:
    """Class for segments of the bandstructure. A Segment contains data for a
    particular region of reciprocal space and particular band.

    Attributes:
        band (int): The band number of the segment (counting starts from 0).
        kpoint_indices (list(int)): The indices of :attr:`effmass.inputs.Data.kpoints` from which this Segment is formed.
        kpoints (array(float)): 2-dimensional array. Each row contains the fractional coordinates of a kpoint [kx,ky,kz]. A slice of :attr:`effmass.inputs.Data.kpoints`.
        cartesian_kpoints (array(float)): 2-dimensional array. Each row contains the cartesian coordinates (angstrom :math:`^{-1}`) of a kpoint.
        dk_angs (array(float)): 1-dimensional array which contains the distance (angstrom :math:`^{-1}`) between each kpoint and the extrema.
        dk_bohr (array(float)): 1-dimensional array which contains the distance (bohr :math:`^{-1}`) between each kpoint and the extrema.
        energies (array(float)): 1-dimensional array which contains the energy (eV) of each eigenstate. A slice of :attr:`effmass.inputs.Data.energies`.
        dE_eV (array(float)): 1-dimensional array which contains the difference in energy (hartree) between each kpoint and the extrema.
        dE_hartree (array(float)): 1-dimensional array which contains the difference in energy (eV) between each kpoint and the extrema.
        occupancy (array(float)): 2-dimensional array. Each row contains occupation number of the eigenstates for a particular band. A slice of :attr:`effmass.inputs.Data.occupancy`.
        direction (array): 1-dimensional array with length 3. The direction between kpoints in the segment.
        band_type (str): The band type, determined by occupancy of the eigenstate. Argument choices are "conduction_band", "valence_band" or "unknown". If set to "unknown", some class methods will raise a warning and return None. 
        fermi_energy (float): the fermi energy in eV.
        dos (array): 2-dimensional array. Each row contains density of states data (units "number of states / unit cell")  at a given energy: [energy(float),dos(float)]. A slice of :attr:`effmass.inputs.Data.dos`.
        integrated_dos (array): 2-dimensional array. Each row contains integrated density of states data at a given energy: [energy(float),integrated_dos(float)]. A slice of :attr:`effmass.inputs.Data.integrated_dos`.
           

    """

    def __init__(self, Data, band, kpoint_indices):
        """Initialise an instance of the Segment class.

        Args:
            Data (Data): Data instance initialised from the bandstructure which contains the segment
            band (int): the band number of the segment
            kpoint_indices (list(int)): the kpoint indices of the segment

        Returns:
            None.
        """

        self.band = band
        self.kpoint_indices = kpoint_indices
        self.kpoints = np.array([Data.kpoints[k] for k in kpoint_indices
                                 ])  # in units 2*pi/angstrom
        self.cartesian_kpoints = np.array([
            np.dot(k, Data.reciprocal_lattice) for k in self.kpoints
        ])  # in units 1/Angstrom. Reciprocal lattice includes factor 2*pi.
        self.dk_angs = np.linalg.norm(
            self.cartesian_kpoints - self.cartesian_kpoints[0], axis=1)
        self.dk_bohr = np.divide(
            self.dk_angs, angstrom_to_bohr
        )  # divide as we are in reciprocal space --> units in inverse length
        self.energies = np.array(
            [Data.energies[band, k] for k in kpoint_indices])  # in units eV
        self.dE_eV = self.energies - self.energies[0]
        self.dE_hartree = np.multiply(self.energies - self.energies[0],
                                      ev_to_hartree)
        self.occupancy = np.array(
            [Data.occupancy[band, k] for k in kpoint_indices])
        self.direction = extrema.calculate_direction(self.kpoints[1],
                                                     self.kpoints[2])
        self.fermi_energy = Data.fermi_energy
        self.dos = self._dos(Data)
        self.integrated_dos = self._integrated_dos(Data)
        self._VBM = Data.VBM
        self._CBM = Data.CBM
        self.band_type = self._band_type()

    def __str__(self):
        """
        Segment string method.

        Returns:
            A string containing the energy of the Segment extrema (referenced to the VBM) and the start- and end- points of the Segment in reciprocal space.
        """
        energy_str = "{0:.2f}".format(self.energies[0] - self._VBM)
        start_str = str(np.round(self.kpoints[0],3))
        end_str = str(np.round(self.kpoints[-1],3))
        return energy_str+" eV; "+start_str+"-->"+end_str

    def _check_kanefit_points(self, polyfit_order=6):
        """Raises an AssertionError if there are not enough data points to fit
        a Kane dispersion.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.

        Returns:
            None.
        """
        idx = self.explosion_index(polyfit_order=polyfit_order)
        assert idx > 3,"Unable to calculate alpha parameter as inflexion too close to extrema"

    def explosion_index(self, polyfit_order=6):
        r"""
        This will find the index at which there is a change in sign of the second derivative. 

        In the region of this point the first derivative will pass through zero and so the transport mass (:math:`\frac{1}{\delta E \delta k}`) will explode.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.

        Notes:
            This marks the point at which the Kane dispersion is definitely not valid, although it may be that the Kane dispersion is a poor approximation prior to this.
        """
        dedk, d2edk2 = self.poly_derivatives(
            polyfit_order=polyfit_order,
            dk=self.dk_bohr,
            polyfit_weighting=False)
        sign = np.sign(d2edk2)
        signchange = ((np.roll(sign, 1) - sign) != 0).astype(int)
        signchange[0] = 0
        if 1 in signchange:
            cutoff = np.where(signchange == 1)[0][0]
        else:
            cutoff = len(self.dE_eV) - 1
        return cutoff

    def _fermi_function(self, eigenvalue, fermi_level=None, temp=300):
        r""" 
        Calculates the probability that an eigenstate is occupied using Fermi-Dirac statistics:
            ..math::
                p=\frac{1}{e^{\frac{Delta E}{kT}+1}

        Args:
            eigenvalue (float): the energy (eV) of the eigenstate.
            fermi_level (float): the fermi level (eV).

        Returns:
            float: the probability that the eigenstate is occupied
        
        """
        if temp <= 0:
            raise ValueError("temperature must be more than 0K")
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level
        if self.band_type == "conduction_band":
            probability = (1 / ((np.exp(
                (eigenvalue - fermi_level) / (((boltzmann * temp) / q))) + 1)))
        elif self.band_type == "valence_band":
            probability = np.subtract(1,(1 / ((np.exp(
                (eigenvalue - fermi_level) / (((boltzmann * temp) / q))) + 1))))
        else:
            raise ValueError("Unable to calculate fermi function as the band type is unknown. Please set the Segment.band_type attribute manually (options are \"valence_band\" or \"conduction_band\").")
        return probability

    def weighting(self, fermi_level=None, temp=300):
        """Calculates a weighting for each kpoint using the Fermi-Dirac
        statistics.

        Args:
            quasi_fermi_level (float, optional): The fermi level to be used in Fermi-Dirac statistics. Defaults to :attr:`effmass.inputs.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.

        Returns:
            array(float): A 1-dimensional array which contains the weighting for each k-point.

        Notes:
            The weighting is relative only: constants will cancel out when using to weight least square fits or means.
        """

        fermi_level = self.fermi_energy if fermi_level is None else fermi_level
        # weighting is given by the fermi function
        weighting = (np.array([
            self._fermi_function(E, fermi_level=fermi_level, temp=temp)
            for E in self.energies
        ]))

        assert weighting.all() != 0,"Unable to assign a weighting to the dispersion as the Fermi-Dirac distribution at all kpoints is smaller than Python's double precision floats. Calculate band edge values using Segment.five_point_leastsq_effmass() or Segment.finite_difference_effmass()."

        return weighting

    def _band_type(self):
        """Identifies the band type, determined by the occupancy of the
        first k-point in the segment.

        An occupancy of 1.0 or 2.0 returns "valence_band" (as :class:`~effmass.analysis.Segment` is in the valence band).
        An occupancy of 0.0 returns "conduction_band" (as :class:`~effmass.analysis.Segment` is in the conduction band).
        In the case of partial occupancy, "unknown" is returned.

        Args:
            None

        Returns:
            str: the band type of the segment, either "valence_band", "conduction_band" or "unknown".
        """

        if self.energies[0] <= self._VBM:
            band_type = "valence_band"
        elif self.energies[0] >= self._CBM:
            band_type = "conduction_band"
        else:
            print(
                "The band type unknown. Please set the Segment.band_type attribute manually."
            )
            band_type = "unknown"
        return band_type

    def _dos(self, Data):
        """Returns slice of :attr:`effmass.Data.dos` corresponding to the
        energy range of the segment.

        Args:
            Data (Data): :class:`~effmass.inputs.Data` instance which was used to generate the :class:`~effmass.analysis.Segment`.

        Returns:
            list(floats): slice of :attr:`effmass.Data.dos` corresponding to the energy range of the segment.
        """
        if Data.dos == []:
            dos = []
        else:
            dos = []
            for i in range(len(self.energies)):
                for j in range(len(Data.dos)):
                    if self.energies[i] < Data.dos[j][0]:
                        dos.append(Data.dos[j][1])
                        break

        return dos

    def _integrated_dos(self, Data):
        """Returns slice of :attr:`effmass.Data.integrated_dos` corresponding
        to the energy range of the segment.

        Args:
            Data (Data): :class:`~effmass.inputs.Data` instance which was used to generate the :class:`~effmass.analysis.Segment`.

        Returns:
            array: slice of :attr:`effmass.Data.integrated_dos` corresponding to the energy range of the segment.
        """
        if Data.integrated_dos == []:
            integrated_dos = []
        else:
            integrated_dos = []
            for i in range(len(self.energies)):
                for j in range(len(Data.integrated_dos)):
                    if self.energies[i] < Data.integrated_dos[j][0]:
                        integrated_dos.append(Data.integrated_dos[j][1])
                        break
        return integrated_dos

    # The collection of methods below calculate the optical effective mass by integrating numerically the analytical
    # expression for the second derivative of the Kane dispersion multiplied by a Fermi-Dirac weighting.

    def mass_integration(self,
                         fermi_level=None,
                         temp=300,
                         alpha=None,
                         mass_bandedge=None,
                         upper_limit=None):
        """Integrates the product of the fermi-dirac distribution, density of
        states and second derivative of kane dispersion along the one-
        dimensional slice of k-space defined by
        :class:`~effmass.analysis.Segment` (up to
        :meth:`~effmass.analysis.Segment.explosion_index`).

        Args:
            fermi_level (float, optional): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            alpha (float, optional): The alpha parameter of the Kane dispersion (hartree$^{-1}$)
            mass_bandedge: The mass at bandedge parameter of the Kane dispersion (units electron mass).
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: The optical effective mass (units of electron mass), defined as the inverse of the second derivative of a kane dispersion, weighted according to occupancy of available eigenstates (the product of density of states and the fermi-dirac distribution).
        Note:
            The sign of the alpha parameter and mass_bandedge are important. If these are negative (as would be expected for the valence band), then they must be passed as negative values to the function.
        """

        alpha = self.alpha() if alpha is None else alpha
        mass_bandedge = self.kane_mass_band_edge(
        ) if mass_bandedge is None else mass_bandedge
        upper_limit = self.dk_bohr[
            self.explosion_index()] if upper_limit is None else upper_limit
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level

        result = integrate.quad(
            self._mass_integrand,
            0,
            upper_limit,
            args=(fermi_level, temp, alpha, mass_bandedge))
        return result

    def _mass_integrand(self, k, fermi_level, temp, alpha, mass_bandedge):
        """Helper function for
        :meth:`~effmass.analysis.Segment.mass_integration`.

        A function for the weighted mass integrand.
        """
        return self.fd(k, fermi_level, temp, alpha, mass_bandedge) * (
            self._second_derivative_kane_dispersion(k, alpha, mass_bandedge))

    def fd(self, k, fermi_level, temp, alpha, mass_bandedge):
        r"""
        Calculates the probability that an eigenstate of momentum k is occupied, using Fermi-Dirac statistics and assuming a Kane dispersion.

        Args:
            k: the momentum of the eigenstate (bohr:math:`^{-1}`).
            fermi_level (float, optional): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            alpha (float, optional): The alpha parameter of the Kane dispersion (hartree$^{-1}$)
            mass_bandedge: The mass at bandedge parameter of the Kane dispersion (units electron mass).

        Returns:
            float: The probability that the eigenstate is occupied.
        Note:
            The sign of the alpha parameter and mass_bandedge are important. If these are negative (as would be expected for the valence band), then they must be passed as negative values to the function.
           
        """
        assert temp > 0, "temperature must be more than 0K"
        if self.band_type == "conduction_band":
            return (1 / ((np.exp((self.energies[0] + self._kane_dispersion(
                k, alpha, mass_bandedge) - fermi_level) / (
                    (boltzmann * temp) / q))) + 1))
        elif self.band_type == "valence_band":
            return 1 - (1 / ((np.exp((self.energies[0] + self._kane_dispersion(
                k, alpha, mass_bandedge) - fermi_level) / (
                    (boltzmann * temp) / q))) + 1))
        elif self.band_type == "unknown":
            raise ValueError("Unable to calculate fermi function as there is partial occupancy of the bands and the band type is unknown. Please set the Segment.band_type attribute manually (options are \"valence_band\" or \"conduction_band\").")
        else:
            raise ValueError("Please set the Segment.band_type attribute (options are \"valence_band\" or \"conduction_band\")")

    def _kane_dispersion(self, k, alpha, mass_bandedge):
        """Helper function for :meth:`~effmass.analysis.Segment.fd`.

        Analytic form of the kane dispersion.
        """
        # returns in eV

        return (_solve_quadratic(
            alpha, 1, [(-k * k) / (2 * mass_bandedge)])[0]) / ev_to_hartree

    def _second_derivative_kane_dispersion(self, k, alpha, mass_bandedge):
        """Helper function for
        :meth:`~effmass.analysis.Segment.mass_integration`.

        Analytic form of the second derivative of the kane dispersion.
        """
        # returns in eV
        return 1 / (mass_bandedge * ((1 + (
            (2 * k * k * alpha) / (mass_bandedge)))**(3 / 2)))

    def weight_integration(self,
                           fermi_level=None,
                           temp=300,
                           alpha=None,
                           mass_bandedge=None,
                           upper_limit=None):
        """Integrates the product of the fermi-dirac distribution
        along the one-dimensional slice of k-space defined by
        :class:`~effmass.analysis.Segment` (up to
        :meth:`~effmass.analysis.Segment.explosion_index`).

        Args:
            fermi_level (float, optional ): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: A normalisation factor for :meth:`~effmass.analysis.Segment.mass_integration`.
        Note:
            The sign of the alpha parameter and mass_bandedge are important. If these are negative (as would be expected for the valence band), then they must be passed as negative values to the function.
       
        """
        alpha = self.alpha() if alpha is None else alpha
        mass_bandedge = self.kane_mass_band_edge(
        ) if mass_bandedge is None else mass_bandedge
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level
        upper_limit = self.dk_bohr[
            self.explosion_index()] if upper_limit is None else upper_limit
        result = integrate.quad(
            self._weight_integrand,
            0,
            upper_limit,
            args=(fermi_level, temp, alpha, mass_bandedge))
        return result

    def _weight_integrand(self, k, fermi_level, temp, alpha, mass_bandedge):
        """Helper function for
        :meth:`~effmass.analysis.Segment.weight_integration`.

        A function for the weight integrand.
        """
        return self.fd(k, fermi_level, temp, alpha, mass_bandedge)

    def optical_effmass_kane_dispersion(self,
                                        fermi_level=None,
                                        temp=300,
                                        alpha=None,
                                        mass_bandedge=None,
                                        upper_limit=None):
        r"""
        Calculates the optical effective mass, with the dispersion approximated by a Kane quasi-linear function.

        This optical effective mass is defined as:
            ..math::
                \frac{1}{m_o} = \frac{\int f(E_k(k),T) \frac{\delta^2 E_k(k)}{\delta k^2} dk}{\int f(E_k(k),T)  dk}

        where the integral is along the one-dimensional slice of k-space defined by :class:`~effmass.analysis.Segment` (up to :meth:`~effmass.analysis.Segment.explosion_index`) and :math:`f(E_k,T)` is the Fermi-Dirac distribution. E_k(k) is the Kane dispersion:
            ..math::
                \frac{\hbar ^2}{2m_{tb}} = E(1+\alpha E)

        where the transport mass at bandedge (:math:`m_{tb}`) is calculated using :meth:`effmass.analysis.Segment.kane_mass_band_edge` and the alpha parameter is calculated using :meth:`effmass.analysis.Segment.alpha`.
        
        Args:
            fermi_level (float, optional ): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            alpha (float, optional): The alpha parameter of the Kane dispersion (hartree$^{-1}$)
            mass_bandedge: The mass at bandedge parameter of the Kane dispersion (units electron mass).
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: The optical effective mass (units of electron mass) of the :class:`~effmass.analysis.Segment`.

        Note:
           The sign of the alpha parameter and mass_bandedge are important. If these are negative (as would be expected for the valence band), then they must be passed as negative values to the function.
        """

        fermi_level = self.fermi_energy if fermi_level is None else fermi_level

        if ((fermi_level) > (self.energies[self.explosion_index()])):
            if (self.band_type == "conduction_band"):
                print(
                    "Fermi level {} is beyond the energy range where the Kane dispersion is valid.".
                    format(fermi_level))

        if (fermi_level) < (self.energies[self.explosion_index()]):
            if self.band_type == "valence_band":
                print(
                    "Fermi level {} is beyond the energy range where the Kane dispersion is valid.".
                    format(fermi_level))

        top = self.mass_integration(
            fermi_level=fermi_level,
            temp=temp,
            alpha=alpha,
            mass_bandedge=mass_bandedge,
            upper_limit=upper_limit)
        bottom = self.weight_integration(
            fermi_level=fermi_level,
            alpha=alpha,
            mass_bandedge=mass_bandedge,
            temp=temp,
            upper_limit=upper_limit)

        assert top[0] != 0, "Unable to calculate the optical effective mass as the Fermi-Dirac distribution at all kpoints is smaller than Python's double precision floats."

        return bottom[0] / top[0]

    ####

    def _polynomial_function(self, polyfit_order=6, polyfit_weighting=True):
        """Helper function which constructs a polynomial function using a least
        squares fit to dispersion data."""
        _check_poly_order(polyfit_order)
        # we know that there is symmetry between k and -k
        negative_dk = [-value for value in self.dk_bohr[::-1]]
        sym_dk = np.concatenate((negative_dk, self.dk_bohr[1:]))
        sym_dE = np.concatenate((self.dE_hartree[::-1], self.dE_hartree[1:]))
        if polyfit_weighting:
            # weight to enable a better fit for the values where it is important
            weighting = self.weighting()
        else:
            weighting = np.ones(len(self.dE_hartree))
        W = np.append(weighting[::-1],
                      weighting[1:])  # as it needs same dimensions as x and y.
        W = np.sqrt(
            np.diag(W))  # to allow dot product between weight and matrix/y.
        # for explanation of the coefficient matrix and least squares fit see https://stackoverflow.com/questions/32260204/numpy-polyfit-passing-through-0
        # for how to incorporate weighting see https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python
        # and https://stackoverflow.com/questions/19624997/understanding-scipys-least-square-function-with-irls
        # by eliminating the units matrix of the array we are forcing a zero offset; the fit must pass through 0,0 as is physical
        coeff_matrix = np.vstack([
            sym_dk**i for i in np.arange(polyfit_order + 1)[polyfit_order:0:-1]
        ]).T
        w_coeff_matrix = np.dot(W, coeff_matrix)
        w_sym_dE = np.dot(sym_dE, W)
        coeff = np.append(np.linalg.lstsq(w_coeff_matrix, w_sym_dE)[0],
                          [0])  # remember to set zeroth power to 0!
        function = np.poly1d(coeff)
        # function = np.poly1d(np.polyfit(sym_dk,sym_dE,polyfit_order,w=W)) ----> simple polyfit call for sanity's sake
        return function

    def poly_derivatives(self,
                         polyfit_order=6,
                         polyfit_weighting=True,
                         dk=None):
        """Constructs a polynomial function using a least squares fit to
        Segment dispersion data, then evaluates first and second order
        derivatives.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate first and second order derivatives. Defaults to 100 points evenly distributed across the whole segment.

        Returns:
            tuple: A tuple containing a 1d array of first derivatives and 1d array of second derivatives, evaluated at dk: ([dedk],[d2edk2])
        """
        function = self._polynomial_function(
            polyfit_order=polyfit_order, polyfit_weighting=polyfit_weighting)
        if dk is None:
            dk = np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)
        dedk = function.deriv(m=1)(dk)
        d2edk2 = function.deriv(m=2)(dk)

        return dedk, d2edk2

    def poly_fit(self, polyfit_order=6, polyfit_weighting=True):
        """Constructs a polynomial function using a least squares fit to
        :class:`~effmass.analysis.Segment` dispersion data, then evaluates at
        100 points evenly distributed across
        :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.

        Returns:
            array: 1d array containing energies (hartree).
        """
        function = self._polynomial_function(
            polyfit_order=polyfit_order, polyfit_weighting=polyfit_weighting)
        values = np.polyval(
            function, np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100))
        return values

    def optical_poly_effmass(self, polyfit_order=6):
        r"""
        Calculates the optical effective mass with a polynomial approximation to the dispersion

        This optical effective mass is defined as:
            ..math::
                \frac{1}{m_o} = \frac{\sum_i f(E_i) g(k_i) \frac{\delta^2 E}{\delta k^2}|_i}{\sum f(E_i) g(k_i)}

        where the sum is over eigenstates i contained withing the :class:`~effmass.analysis.Segment`. :math:`f_(E_i)` is the probability of occupation (Fermi-Dirac statistics) and :math:`g(k_i)` is the density of states at that k-point.
        
        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
        
        Returns:
            float: The optical effective mass (units of electron mass) of the :class:`~effmass.analysis.Segment`.
        """
        inertial_mass = self.inertial_effmass(
            polyfit_order=polyfit_order, dk=self.dk_bohr)
        weighting = self.weighting()
        optical_mass = np.power(
            np.average(np.power(inertial_mass, -1), weights=weighting), -1)
        return optical_mass

    def inertial_effmass(self,
                         polyfit_order=6,
                         dk=None,
                         polyfit_weighting=False):
        r"""
        Calculates the inertial (curvature) effective mass (:math:`\frac{1}{\frac{\delta^2E}{\delta k^2}}`), evaluated at :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate second order derivatives. Defaults to 100 points evenly distributed across the whole segment.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
        
        Returns:
            array(float). 1d array containing the conductivity effective mass (in units of electron rest mass) evaluated at the points specified in dk.
        """
        dedk, d2edk2 = self.poly_derivatives(
            polyfit_order=polyfit_order,
            dk=dk,
            polyfit_weighting=polyfit_weighting)
        conductivity_mass = [(1 / x) for x in d2edk2]
        return conductivity_mass

    def transport_effmass(self,
                          polyfit_order=6,
                          dk=None,
                          polyfit_weighting=False):
        r"""
        Calculates the transport mass (:math:`\frac{k}{\delta E \delta k}`), evaluated at :attr:`~effmass.analysis.Segment.dk_bohr` .

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate first order derivatives. Defaults to 100 points evenly distributed across the whole segment.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
           
        Returns:
            array(float). 1d array containing the transport effective mass (in units of electron rest mass) evaluated at the points specified in dk.
        """
        dedk, d2edk2 = self.poly_derivatives(
            polyfit_order=polyfit_order,
            dk=dk,
            polyfit_weighting=polyfit_weighting)
        # dk=0 for first term gives discontinuity
        if dk is None:
            dk = np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)
        transport_mass = [(x / y) for x, y in zip(dk[1:], dedk[1:])]
        transport_mass = np.append(np.array([np.nan]), transport_mass)
        return transport_mass

    def alpha(self, polyfit_order=6, truncate=True):
        r"""
        The transport mass (:math:`\frac{k}{\delta E \delta k}`) is calculated as a function of :attr:`~effmass.analysis.Segment.dk_bohr` and fitted to a straight line. The gradient of this line determines the alpha parameter which is used in the kane dispersion.
        
        See :meth:`effmass.analysis.Segment.transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            truncate (bool, optional): If True, data only up to :meth:`effmass.analysis.Segment.explosion_index` is used. If False, alpha is calculated using data for the whole segment. Defaults to True.

        Returns:
            float: The alpha parameter (hartree :math:`^{-1}`).
        """
        if truncate:
            self._check_kanefit_points(polyfit_order=polyfit_order)
        transport_mass = self.transport_effmass(
            polyfit_order=polyfit_order,
            dk=self.dk_bohr,
            polyfit_weighting=False)
        if truncate:
            idx = self.explosion_index(polyfit_order=polyfit_order)
            gradient, intercept = np.polyfit(self.dE_hartree[1:idx + 1],
                                             transport_mass[1:idx + 1], 1)
        else:
            gradient, intercept = np.polyfit(self.dE_hartree[1:],
                                             transport_mass[1:], 1)
        alpha = np.divide(gradient, 2 * intercept)
        return alpha

    def kane_mass_band_edge(self, polyfit_order=6, truncate=True):
        r"""
        The transport mass (:math:`\frac{1}{\delta E \delta k}`) is calculated as a function of :attr:`~effmass.analysis.Segment.dk_bohr` and fitted to a straight line. The intercept of this line with the y-axis gives a transport mass at bandedge which is used as a parameter in the kane dispersion.
        
        See :meth:`effmass.analysis.Segment.transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            truncate (bool, optional): If True, data only up to :meth:`effmass.analysis.Segment.explosion_index` is used. If False, alpha is calculated using data for the whole segment. Defaults to True.
        
        Returns:
            float: transport mass at bandedge (in units of electron mass).
        """
        if truncate:
            self._check_kanefit_points(polyfit_order=polyfit_order)
        transport_mass = self.transport_effmass(
            polyfit_order=polyfit_order,
            dk=self.dk_bohr,
            polyfit_weighting=False)
        if truncate:
            idx = self.explosion_index(polyfit_order=polyfit_order)
            gradient, intercept = np.polyfit(self.dE_hartree[1:idx + 1],
                                             transport_mass[1:idx + 1], 1)
        else:
            gradient, intercept = np.polyfit(self.dE_hartree[1:],
                                             transport_mass[1:], 1)
        return intercept  # intercept is the band edge transport mass

    def kane_fit(self, polyfit_order=6):
        r"""
        Calculate the Kane quasi-linear dispersion parameters, then evaluates at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.

        The Kane quasi-linear dispersion is described by 
            ..math::
                \frac{\hbar ^2 k^2}{2m_{tb}} = E(1+\alpha E)

        where the transport mass at bandedge (:math:`m_{tb}`) and the alpha parameter are calculated by fitting a linear function to the transport mass :math:`m_t` as a function of energy E.
            ..math::
                m_t = m_{tb}(1+2\alpha E)

        The transport mass :math:`m_t` is calculated by approximating the dispersion with a polynomial function and taking the first derivative, see :meth:`~effmass.analysis.Segment.transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
 
        Returns:
            array: 1d array containing energies (hartree).
        """
        self._check_kanefit_points(polyfit_order=polyfit_order)
        alpha = self.alpha(polyfit_order=polyfit_order)
        transport_mass_bandedge = self.kane_mass_band_edge(
            polyfit_order=polyfit_order)
        bandedge_energy = [
            ((k**2)) / (2 * transport_mass_bandedge)
            for k in np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)
        ]
        roots = _solve_quadratic(alpha, 1,
                                 [(-1) * ele for ele in bandedge_energy])
        return roots

    def finite_difference_effmass(self):
        """The curvature at the band edge is calculated using a second order
        forward finite difference method. This is then inverted to give an
        effective mass.

        Args:
            None

        Returns:
            float: Bandedge effective mass from finite difference (in units of electron mass).
        """
        x = (self.dE_hartree[2] -
             (2 * self.dE_hartree[1])) / (self.dk_bohr[1]**2)
        mass = 1 / x
        return mass

    def finite_difference_fit(self):
        r"""
        Calculates the curvature at the band edge using a finite difference method and then evaluates the corresponding quadratic dispersion along :attr:`~effmass.analysis.Segment.dk_bohr`.

        See :meth:`effmass.analysis.Segment.finite_difference_effmass`.

        Args:
            None

        Returns:
            list(float): list containing energies (hartree). The energies are calculated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr` using the quadratic approximation. 
        """
        m_bandedge = self.finite_difference_effmass()
        values = [((k**2) / (2 * m_bandedge))
                  for k in np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)]
        return values

    def weighted_leastsq_effmass(self):
        """Fits a parabolic dispersion using the weighted least-squares method
        to all points in the segment, plus those from symmetry, E(k)=E(-k).

        Args:
            None

        Returns:
            float: Curvature effective mass (in units of electron mass)

        Notes:
            weighting is given by the Fermi-Dirac distribution.
        """

        # https://stackoverflow.com/questions/19624997/understanding-scipys-least-square-function-with-irls
        # https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python

        negative_dk = [-value for value in self.dk_bohr[::-1]]
        sym_dk = np.concatenate((negative_dk, self.dk_bohr[1:]))
        coeff_matrix = np.vstack([sym_dk**2]).T

        weighting = self.weighting()
        W = np.append(weighting[::-1],
                      weighting[1:])  # as it needs same dimensions as x and y.
        W = np.sqrt(
            np.diag(W))  # to allow dot product between weight and matrix/y.

        sym_dE = np.concatenate((self.dE_hartree[::-1], self.dE_hartree[1:]))
        w_sym_dE = np.dot(sym_dE, W)
        w_coeff_matrix = np.dot(W, coeff_matrix)

        coeff = np.linalg.lstsq(w_coeff_matrix, w_sym_dE)[0]
        mass = 1 / (2 * coeff[0])
        return mass

    def weighted_leastsq_fit(self):
        """Calculates the curvature effective mass using a weighted least-
        squares fit and then evaluates the corresponding parabolic dispersion
        along :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            None

        Returns:
            list(float): list containing energies (hartree). The energies are calculated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.
        """

        m_bandedge = self.weighted_leastsq_effmass()
        values = [((k**2) / (2 * m_bandedge))
                  for k in np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)]
        return values

    def five_point_leastsq_effmass(self):
        """Fits a parabolic dispersion using the least-squares method to 5
        points (3 DFT-calculated points + 2 from symmetry).

        Args:
            None

        Returns:
            float: Curvature effective mass (in units of electron mass).

        Notes:
            no weighting is used.
        """

        sym_dk = np.array([
            -self.dk_bohr[2], -self.dk_bohr[1], self.dk_bohr[0],
            self.dk_bohr[1], self.dk_bohr[2]
        ])
        sym_dE = np.array([
            self.dE_hartree[2], self.dE_hartree[1], self.dE_hartree[0],
            self.dE_hartree[1], self.dE_hartree[2]
        ])
        coeff_matrix = np.vstack([sym_dk**2]).T
        coeff = np.linalg.lstsq(coeff_matrix, sym_dE)[0]
        mass = 1 / (2 * coeff[0])
        return mass

    def five_point_leastsq_fit(self):
        """Calculates the curvature effective mass using a parabolic least-
        squares fit and then evaluates the corresponding parabolic dispersion
        along :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            None

        Returns:
            list(float): list containing energies (hartree). The energies are calculated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.
        """
        m_bandedge = self.five_point_leastsq_effmass()
        values = [((k**2) / (2 * m_bandedge))
                  for k in np.linspace(self.dk_bohr[0], self.dk_bohr[-1], 100)]
        return values


