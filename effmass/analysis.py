#! /usr/bin/env python3

import warnings
import itertools
import scipy
import numpy as np
from scipy import integrate
from effmass import extrema
from effmass import angstrom_to_bohr
from effmass import ev_to_hartree
from effmass import boltzmann
from effmass import q

warnings.simplefilter(action='ignore', category=FutureWarning)

def _check_poly_order(polyfit_order):
    """ Raises an AssertionError if the order of the polynomial is less than 2.

    Args:
        polyfit_order (int): the order of the polynomial fit.

    Returns:
        None.
    """
    assert polyfit_order > 1, "the order of the polynomial must be 2 or higher"

def weighted_mean(values, weighting):
    r"""
    Calculates the mean of values :math:`v_i` with weighting :math:`w_i`: :math:`\frac{\sum_i(v_iw_i)}{\sum_i(w_i)}`.

    Args:
        values (list(float)): values of which we want the weighted mean.
        weighting (list(float)): weighting for each of these values.
    
    Returns:
        float: the weighted mean of the values.
    
    """
    weighted_mean= np.sum(np.multiply(values, weighting)) / np.sum(weighting)
    return weighted_mean

def solve_quadratic(a,b,c):
    r"""
    Solves quadratic equation of the form :math:`ax^2+bx+c=0` for each value of c given.

    Args:
        a (float): the coefficient of :math:`x^2`.
        b (float): the coefficient of :math:`x^1`.
        c (list(float)): a list of coefficients of :math:`x^0`.

    Returns:
        list(float): the root of the equation for each value of c.
    
    """
    det = [b**2 - 4*a*item for item in c]
    x = []
    for d in det:
        if d < 0:
            x.append(np.nan)
        if d == 0:
            x.append(-b/(2*a))
        if d > 0:
            x.append((-b+np.sqrt(d))/(2*a))
    return x

def _average_allband_optical_effmass(group,polyfit_order=6, optical_weighted_average=True):
    """
    Helper function for calc_allband_optical_effmass.

    Args:
        group (iterator): iterator of :class:`~effmass.analysis.Segment` instances grouped by direction and occupation.
        polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
        optical_weighted_average (bool, optional): If True, a weighting (proportional to the product of the Fermi-Dirac function and density of states) is applied. If False, no weighting is applied. Defaults to True. 
        use_doscar (bool, optional): If True, DOSCAR data is used as density of states weighting. If False, :math:`k^2` is used. Defaults to False.

    Returns:
        float: optical mass calculated using multiple segments which have the same direction and occupation. 
    """

    mass = []
    weight = []
    for segment in group:
        mass.append(segment.calc_conductivity_effmass(polyfit_order=polyfit_order, dk=segment.dk_bohr))
        if optical_weighted_average:
            weight.append(segment.return_weighting())
        else:
            weight.append(np.ones(len(segment.energies)))
    optical_mass = np.power(weighted_mean(np.power(mass[0],-1), weight[0]),-1)
    return optical_mass

def calc_allband_optical_effmass(segments,polyfit_order=6,optical_weighted_average=True):
    r"""
    Calculates the optical effective mass across multiple segments, where the dispersion of each Segment has been approximated with a polynomial function.

    This optical effective mass :math:`m_o` is defined as:
        .. math::
           \frac{1}{m_o} = \frac{ \sum_b \sum_i (f(E_i)g(k_i)\frac{\delta^2E}{\delta k^2}|_i)}{\sum_b \sum_i (f(E_i)g(k))}

    where :math:`\frac{\delta^2E}{\delta k^2}` is the second derivative of the polynomial fit evaluated at :attr:`~effmass.analysis.Segment.dk_bohr`, :math:`\sum_b` is the sum over all bands which are in the same direction and have the same occupancy,
    :math:`\sum_i` is the sum over all eigenstates within each band Segment, :math:`f_(E_i)` is the probability of occupation (Fermi-Dirac statistics),
    :math:`g(k_i)` is the density of states at that k-point.
    Segments can be in different directions and have different occupancy.
    The function will sort them accordingly and average across only those which have common direction and occupation.

    Args:
        segments (list(Segment)): A list of Segment instances.
        polyfit_order (int, optional): order of polynomial used to approximate the dispersion.
            Defaults to 6.
        optical_weighted_average (bool, optional): If True, a weighting (proportional to the product of the Fermi-Dirac function and density of states) is applied. If False, no weighting is applied. Defaults to True. See :meth:`effmass.analysis.Segment.return_weighting`.
        use_doscar (bool, optional): If True, :attr:`effmass.inputs.Data.dos` is used as density of states weighting. If False, :attr:`~effmass.analysis.Segment.dk_bohr`:math:`^2` is used. Defaults to False. See :meth:`effmass.analysis.Segment.return_weighting`.

    Returns:
        list: A 2-dimensional list. Each row contains [:attr:`~effmass.analysis.Segment.direction`,:attr:`~effmass.analysis.Segment.occupancy`,optical_mass,iterator] for each group of Segment instances with share the same occupation and direction.
        The optical mass is an average across the group of Segment instances and is given in units of electron mass. The iterator is an <itertools group>`https://docs.python.org/2/library/itertools.html#itertools.groupby`_ object.

    """
    allband_effmass = []
    sorted_segments = sorted(segments, key=lambda x: str(np.absolute(x.direction)))
    for direction, group in itertools.groupby(sorted_segments,lambda x: str(np.absolute(x.direction))):
        sorted_group = sorted(group, key=lambda x: x.occupancy[0])
        for occupation, group in itertools.groupby(sorted_group,lambda x: x.occupancy[0]):    
            optical_mass = _average_allband_optical_effmass(group,polyfit_order=polyfit_order,optical_weighted_average=optical_weighted_average)
            allband_effmass.append([direction, occupation, optical_mass, group])

    return allband_effmass


class Segment:
    """
    Class for segments of the bandstructure. A Segment contains data for a particular region of reciprocal space and particular band.

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
        ptype (str): The quasi particle type, determined by occupancy of the eigenstate.
        bandedge_energy: The enery of the VBM (if Segment instance is in valence band) or CBM (if Segment instance is in conduction band). Units are eV.
        fermi_energy (float): the fermi energy in eV. 
        dos (array): 2-dimensional array. Each row contains density of states data (units "number of states / unit cell")  at a given energy: [energy(float),dos(float)]. A slice of :attr:`effmass.inputs.Data.dos`.
        integrated_dos (array): 2-dimensional array. Each row contains integrated density of states data at a given energy: [energy(float),integrated_dos(float)]. A slice of :attr:`effmass.inputs.Data.integrated_dos`.
    """

    def __init__(self, Data, band, kpoint_indices):
        """
        Initialise an instance of the Segment class.

        Args:
            Data (Data): Data instance initialised from the bandstructure which contains the segment
            band (int): the band number of the segment
            kpoint_indices (list(int)): the kpoint indices of the segment

        Returns:
            None.
        """


        self.band = band
        self.kpoint_indices = kpoint_indices
        self.kpoints = np.array([Data.kpoints[k] for k in kpoint_indices]) # in units 2*pi/angstrom 
        self.cartesian_kpoints = np.array([np.dot(k, Data.reciprocal_lattice) for k in self.kpoints]) # in units 1/Angstrom. Reciprocal lattice includes factor 2*pi.
        self.dk_angs = np.linalg.norm(self.cartesian_kpoints - self.cartesian_kpoints[0], axis=1)
        self.dk_bohr = np.divide(self.dk_angs,angstrom_to_bohr) # divide as we are in reciprocal space --> units in inverse length
        self.energies = np.array([Data.energies[band,k] for k in kpoint_indices]) # in units eV
        self.dE_eV = self.energies - self.energies[0]
        self.dE_hartree = np.multiply(self.energies - self.energies[0],ev_to_hartree)
        self.occupancy = np.array([Data.occupancy[band,k] for k in kpoint_indices])
        self.direction = extrema.calculate_direction(self.kpoints[1],self.kpoints[2])
        self.ptype = self.return_ptype()
        self.bandedge_energy = self.return_bandedge_energy(Data)
        self.fermi_energy = Data.fermi_energy
        self.dos = self.return_dos(Data)
        self.integrated_dos = self.return_integrated_dos(Data)


    def fermi_function(self,eigenvalue,fermi_level=None,temp=300):
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
        assert temp>0, "temperature must be more than 0K"
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level
        if self.ptype == "electron":
            probability =  (1 / ((np.exp((eigenvalue-fermi_level)/(((boltzmann*temp)/q))) + 1))) 
        if self.ptype == "hole":
            probability =  1 - (1 / ((np.exp((eigenvalue-fermi_level)/(((boltzmann*temp)/q))) + 1)))
        if self.ptype == "unknown":
            print ("unable to calculate fermi function as particle type unknown (partial occupancy, please set Segment.type)")
            probability = None

        return probability

    def return_ptype(self):
        """
        Identifies the quasi particle, determined by the occupancy of the first k-point in the segment.

        An occupancy of 1.0 or 2.0 returns "hole" (as :class:`~effmass.analysis.Segment` is in the valence band).
        An occupancy of 0.0 returns "electron" (as :class:`~effmass.analysis.Segment` is in the conduction band).
        In the case of partial occupancy, "unknown" is returned.
        
        Args:
            None

        Returns:
            str: the particle type of the segment, either "hole", "electron" or "unknown".
        """

        if self.occupancy[0] == 1.0 or self.occupancy[0]==2.0:
            particle = "hole"
        elif self.occupancy[0] == 0.0:
            particle = "electron"
        else:
            print ("partial occupancy, particle type unknown. Please set Segment.ptype manually.")
            particle = "unknown"
        return particle

    def return_bandedge_energy(self,Data):
        """
        Identifies the bandedge energy of the :class:`~effmass.analysis.Segment`, determined by the occupancy of the first k-point.

        If :attr:`~effmass.analysis.Segment.ptype` is "hole", function returns the Data.VBM (as :class:`~effmass.analysis.Segment` is in the valence band).
        an :attr:`~effmass.analysis.Segment.ptype` is "electron", function returns Data.CBM (as :class:`~effmass.analysis.Segment` is in the conduction band).

        Args:
            Data (Data): :class:`~effmass.inputs.Data` instance which was used to generate the :class:`~effmass.analysis.Segment`.

        Returns:
            float: the bandedge energy of the segment (eV).
        """

        if self.ptype == "hole":
            bandedge_energy =  Data.VBM
        elif self.ptype == "electron":
            bandedge_energy = Data.CBM
        elif self.ptype == "unknown":
            print ("cannot determine bandedge energy as particle type unknown. Please set Segment.ptype.")
            bandedge_energy = None

        return bandedge_energy

    def return_dos(self,Data):
        """
        Returns slice of :attr:`effmass.Data.dos` corresponding to the energy range of the segment.

        Args:
            Data (Data): :class:`~effmass.inputs.Data` instance which was used to generate the :class:`~effmass.analysis.Segment`.

        Returns: 
            list(floats): slice of :attr:`effmass.Data.dos` corresponding to the energy range of the segment.
        """
        if Data.dos == []:
            dos = []
        else:
            dos=[]
            for i in range(len(self.energies)):
                for j in range(len(Data.dos)):
                    if self.energies[i] < Data.dos[j][0]:
                        dos.append(Data.dos[j][1])
                        break

        return dos

    def return_integrated_dos(self,Data):
        """
        Returns slice of :attr:`effmass.Data.integrated_dos` corresponding to the energy range of the segment.

        Args:
            Data (Data): :class:`~effmass.inputs.Data` instance which was used to generate the :class:`~effmass.analysis.Segment`.

        Returns: 
            array: slice of :attr:`effmass.Data.integrated_dos` corresponding to the energy range of the segment.
        """
        if Data.integrated_dos == []:
            integrated_dos = []
        else:
            integrated_dos=[]
            for i in range(len(self.energies)):
                for j in range(len(Data.integrated_dos)):
                    if self.energies[i] < Data.integrated_dos[j][0]:
                        integrated_dos.append(Data.integrated_dos[j][1])
                        break
        return integrated_dos
        
    def return_weighting(self, fermi_level=None, temp=300):
        """
        Calculates a weighting for each kpoint using the Fermi-Dirac statistics.

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
        weighting = (np.array([self.fermi_function(E, fermi_level=fermi_level,temp=temp) for E in self.energies]))  

        return weighting

     #### The collection of methods below calculate the optical effective mass by integrating numerically the analytical
     #### expression for FD,DOS and second derivative of dispersion (currently only the Kane dispersion is implemented).

    def return_mass_integration(self,fermi_level=None,temp=300,alpha=None,mass_bandedge=None,upper_limit=None):
        """
        Integrates the product of the fermi-dirac distribution, density of states and second derivative of kane dispersion along the one-dimensional slice of k-space defined by :class:`~effmass.analysis.Segment` (up to :meth:`~effmass.analysis.Segment.explosion_index`). 

        Args:
            fermi_level (float, optional): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            alpha (float, optional): The alpha parameter of the Kane dispersion (hartree$^{-1}$)
            mass_bandedge: The mass at bandedge parameter of the Kane dispersion (units electron mass).
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: The optical effective mass (units of electron mass), defined as the inverse of the second derivative of a kane dispersion, weighted according to occupancy of available eigenstates (the product of density of states and the fermi-dirac distribution).
        """
        
        alpha = self.calc_alpha() if alpha is None else alpha
        mass_bandedge = self.calc_kane_mass_bandedge() if mass_bandedge is None else mass_bandedge
        upper_limit = self.dk_bohr[self.explosion_index()] if upper_limit is None else upper_limit
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level

        
        result = scipy.integrate.quad(self._mass_integrand,0,upper_limit,args=(fermi_level,temp,alpha,mass_bandedge))
        return result

    def _mass_integrand(self,k,fermi_level, temp,alpha,mass_bandedge):
        """
        Helper function for :meth:`~effmass.analysis.Segment.return_mass_integration`. A function for the weighted mass integrand.
        """
        return self.fd(k,fermi_level,temp,alpha,mass_bandedge)*(self._second_derivative_kane_dispersion(k,alpha,mass_bandedge))

    def fd(self,k,fermi_level,temp,alpha,mass_bandedge):
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
        """
        assert temp>0, "temperature must be more than 0K"
        if self.ptype == "electron":
            return (1 / ((np.exp((self.energies[0]+self._kane_dispersion(k,alpha,mass_bandedge)-fermi_level)/((boltzmann*temp)/q))) + 1)) 
        if self.ptype == "hole":
            return 1 - (1 / ((np.exp((self.energies[0]+self._kane_dispersion(k,alpha,mass_bandedge)-fermi_level)/((boltzmann*temp)/q))) + 1))
        if self.ptype == "unknown":
            print ("unable to calculate fermi function as particle type unknown (partial occupancy, please set Segment.ptype)")
            return
        
    def _kane_dispersion(self,k,alpha,mass_bandedge):
        """
        Helper function for :meth:`~effmass.analysis.Segment.fd`. Analytic form of the kane dispersion.
        """
        # returns in eV

        return (solve_quadratic(alpha,1,[(-k*k)/(2*mass_bandedge)])[0])/ev_to_hartree

    def _second_derivative_kane_dispersion(self,k,alpha,mass_bandedge):
        """
        Helper function for :meth:`~effmass.analysis.Segment.return_mass_integration`. Analytic form of the second derivative of the kane dispersion.
        """
        # returns in eV
        return 1 / (mass_bandedge*((1+ ((2*k*k*alpha)/(mass_bandedge)))**(3/2)))

    def return_weight_integration(self,fermi_level=None,temp=300,alpha=None,mass_bandedge=None,upper_limit=None):
        """
        Integrates the product of the fermi-dirac distribution and density of states along the one-dimensional slice of k-space defined by :class:`~effmass.analysis.Segment` (up to :meth:`~effmass.analysis.Segment.explosion_index`). 

        Args:
            fermi_level (float, optional ): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: A normalisation factor for :meth:`~effmass.analysis.Segment.return_mass_integration`.
        """
        alpha = self.calc_alpha() if alpha is None else alpha
        mass_bandedge = self.calc_kane_mass_bandedge() if mass_bandedge is None else mass_bandedge
        fermi_level = self.fermi_energy if fermi_level is None else fermi_level
        upper_limit = self.dk_bohr[self.explosion_index()] if upper_limit is None else upper_limit
        result = scipy.integrate.quad(self._weight_integrand,0,upper_limit,args=(fermi_level,temp,alpha,mass_bandedge))
        return result

    def _weight_integrand(self,k,fermi_level,temp,alpha,mass_bandedge):
        """
        Helper function for :meth:`~effmass.analysis.Segment.return_weight_integration`. A function for the weight integrand.
        """
        return self.fd(k,fermi_level,temp,alpha,mass_bandedge)

    def calc_optical_effmass_kane_dispersion(self,fermi_level=None,temp=300,alpha=None,mass_bandedge=None,upper_limit=None):
        r"""
        Calculates an optical effective mass for a single :class:`~effmass.analysis.Segment` where the dispersion is approximated with a quasi-linear function.

        This optical effective mass is defined as:
            ..math::
                \frac{1}{m_o} = \frac{\int f(E_k(k),T) \frac{\delta^2 E_k(k)}{\delta k^2} dk}{\int f(E_k(k),T)  dk}

        where the integral is along the one-dimensional slice of k-space defined by :class:`~effmass.analysis.Segment` (up to :meth:`~effmass.analysis.Segment.explosion_index`) and :math:`f(E_k,T)` is the Fermi-Dirac distribution. E_k(k) is the Kane dispersion:
            ..math::
                \frac{\hbar ^2}{2m_{tb}} = E(1+\alpha E)

        where the transport mass at bandedge (:math:`m_{tb}`) is calculated using :meth:`effmass.analysis.Segment.calc_kane_mass_bandedge` and the alpha parameter is calculated using :meth:`effmass.analysis.Segment.calc_alpha`.
        
        Args:
            fermi_level (float, optional ): Fermi level (eV) to be used in Fermi-dirac statistics. Defaults to :attr:`~effmass.analysis.Segment.fermi_energy`.
            temp (float, optional): The temperature (K) to be used in Fermi-Dirac statistics. Defaults to 300.
            alpha (float, optional): The alpha parameter of the Kane dispersion (hartree$^{-1}$)
            mass_bandedge: The mass at bandedge parameter of the Kane dispersion (units electron mass).
            upper_limit (float, optional): The integration upper limit (bohr$^{-1}$). Defaults to where the Kane quasi-linear dispersion is no longer valide, defined by :meth:`~effmass.analysis.Segment.explosion_index`.

        Returns:
            float: The optical effective mass (units of electron mass) of the :class:`~effmass.analysis.Segment`.
        """

        fermi_level = self.fermi_energy if fermi_level is None else fermi_level

        if ((fermi_level) > (self.energies[self.explosion_index()])):
            if (self.ptype == "electron"):
              print ("Fermi level {} is beyond the energy range where the Kane dispersion is valid.".format(fermi_level))


        if (fermi_level) < (self.energies[self.explosion_index()]):
            if self.ptype == "hole":
              print ("Fermi level {} is beyond the energy range where the Kane dispersion is valid.".format(fermi_level))

        top = self.return_mass_integration(fermi_level=fermi_level,temp=temp,alpha=alpha,mass_bandedge=mass_bandedge,upper_limit=upper_limit)
        bottom = self.return_weight_integration(fermi_level=fermi_level,alpha=alpha,mass_bandedge=mass_bandedge,temp=temp,upper_limit=upper_limit)

        return bottom[0]/top[0]

    ####

    def _construct_polynomial_function(self, polyfit_order=6, polyfit_weighting=True):
        """
        Helper function which constructs a polynomial function using a least squares fit to dispersion data.
        """
        _check_poly_order(polyfit_order)
        # we know that there is symmetry between k and -k
        negative_dk = [-value for value in self.dk_bohr[::-1]]
        sym_dk = np.concatenate((negative_dk,self.dk_bohr[1:]))
        sym_dE = np.concatenate((self.dE_hartree[::-1],self.dE_hartree[1:]))
        if polyfit_weighting:
            # weight to enable a better fit for the values where it is important
            weighting=self.return_weighting()
        else:
            weighting = np.ones(len(self.dE_hartree))
        W = np.append(weighting[::-1],weighting[1:]) # as it needs same dimensions as x and y.
        W = np.sqrt(np.diag(W)) # to allow dot product between weight and matrix/y.
        # for explanation of the coefficient matrix and least squares fit see https://stackoverflow.com/questions/32260204/numpy-polyfit-passing-through-0
        # for how to incorporate weighting see https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python
        # and https://stackoverflow.com/questions/19624997/understanding-scipys-least-square-function-with-irls
        # by eliminating the units matrix of the array we are forcing a zero offset; the fit must pass through 0,0 as is physical
        coeff_matrix = np.vstack([sym_dk ** i for i in np.arange(polyfit_order+1)[polyfit_order:0:-1]]).T 
        w_coeff_matrix = np.dot(W,coeff_matrix)
        w_sym_dE = np.dot(sym_dE, W)
        coeff = np.append(np.linalg.lstsq(w_coeff_matrix, w_sym_dE)[0],[0]) # remember to set zeroth power to 0!
        function = np.poly1d(coeff)
        #function = np.poly1d(np.polyfit(sym_dk,sym_dE,polyfit_order,w=W)) ----> simple polyfit call for sanity's sake
        return function
    
    def calc_poly_derivatives(self, polyfit_order=6,polyfit_weighting=True,dk=[]):
        """
        Constructs a polynomial function using a least squares fit to Segment dispersion data, then evaluates first and second order derivatives.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate first and second order derivatives. Defaults to 100 points evenly distributed across the whole segment.
        
        Returns:
            tuple: A tuple containing a 1d array of first derivatives and 1d array of second derivatives, evaluated at dk: ([dedk],[d2edk2])
        """
        function = self._construct_polynomial_function(polyfit_order=polyfit_order, polyfit_weighting=polyfit_weighting)
        if dk == []:
            dk = np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)
        dedk = function.deriv(m=1)(dk)
        d2edk2 = function.deriv(m=2)(dk)

        return dedk, d2edk2

    def return_polyfit(self, polyfit_order=6,polyfit_weighting=True):
        """
        Constructs a polynomial function using a least squares fit to :class:`~effmass.analysis.Segment` dispersion data, then evaluates at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
        
        Returns:
            array: 1d array containing energies (hartree).
        """
        function = self._construct_polynomial_function(polyfit_order=polyfit_order, polyfit_weighting=polyfit_weighting)
        values = np.polyval(function,np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100))
        return values

    def check_kanefit_points(self,polyfit_order=6):
        """
        Raises an AssertionError if there are not enough data points to fit a Kane dispersion. 

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
        
        Returns:
            None.
        """
        idx = self.explosion_index(polyfit_order=polyfit_order)
        assert idx>3, "Unable to calculate alpha parameter as inflexion too close to extrema"

    def return_kanefit(self, polyfit_order=6):
        r"""
        Constructs a kane quasi linear dispersion, then evaluates at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.

        The Kane quasi-linear dispersion is described by 
            ..math::
                \frac{\hbar ^2 k^2}{2m_{tb}} = E(1+\alpha E)

        where the transport mass at bandedge (:math:`m_{tb}`) and the alpha parameter are calculated by fitting a linear function to the transport mass :math:`m_t` as a function of energy E.
            ..math::
                m_t = m_{tb}(1+2\alpha E)

        The transport mass :math:`m_t` is calculated by approximating the dispersion with a polynomial function and taking the first derivative, see :meth:`~effmass.analysis.Segment.calc_transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
 
        Returns:
            array: 1d array containing energies (hartree).
        """
        self.check_kanefit_points(polyfit_order=polyfit_order)
        alpha = self.calc_alpha(polyfit_order=polyfit_order)
        transport_mass_bandedge = self.calc_kane_mass_bandedge(polyfit_order=polyfit_order)
        bandedge_energy = [((k**2)) / (2*transport_mass_bandedge) for k in np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)]
        roots = solve_quadratic(alpha, 1, [(-1)*ele for ele in bandedge_energy]) 
        return roots    

    def return_quadfit(self, polyfit_order=6,polyfit_weighting=True):
        """
        Calculates the curvature at bandedge uing a least squares fit to :class:`~effmass.analysis.Segment` dispersion data. The corresponding pure quadratic function is then evaluated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used. Defaults to True. See :meth:`effmass.analysis.Segment.return_weighting`.
        
        Returns:
            array: 1d array containing energies (hartree).
        """
        dedk, d2edk2 = self.calc_poly_derivatives(polyfit_order=polyfit_order,polyfit_weighting=polyfit_weighting)
        m_bandedge =  1 / d2edk2[0]
        values = [((k**2)/(2*m_bandedge)) for k in np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)]
        return values

    def calc_optical_effmass(self, polyfit_order=6, optical_weighted_average=True):
        r"""
        Calculates an optical effective mass for a single :class:`~effmass.analysis.Segment` where the dispersion is approximated with a polynomial.

        This optical effective mass is defined as:
            ..math::
                \frac{1}{m_o} = \frac{\sum_i f(E_i) g(k_i) \frac{\delta^2 E}{\delta k^2}|_i}{\sum f(E_i) g(k_i)}

        where the sum is over eigenstates i contained withing the :class:`~effmass.analysis.Segment`. :math:`f_(E_i)` is the probability of occupation (Fermi-Dirac statistics) and :math:`g(k_i)` is the density of states at that k-point.
        
        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            optical_weighted_average (bool, optional): If True, a weighting (proportional to the product of the Fermi-Dirac function and density of states) is applied. If False, no weighting is applied. Defaults to True. See :meth:`effmass.analysis.Segment.return_weighting`.
        
        Returns:
            float: The optical effective mass (units of electron mass) of the :class:`~effmass.analysis.Segment`.
        """
        conductivity_mass = self.calc_conductivity_effmass(polyfit_order=polyfit_order,dk=self.dk_bohr)
        if optical_weighted_average:
            weighting = self.return_weighting()
        else:
            weighting = np.ones(len(self.energies))
        optical_mass = np.power(weighted_mean(np.power(conductivity_mass,-1), weighting),-1)
        return optical_mass

    def calc_conductivity_effmass(self, polyfit_order=6,dk=[],polyfit_weighting=False):
        r"""
        Calculates the conductivity mass (:math:`\frac{1}{\frac{\delta^2E}{\delta k^2}}`), evaluated at :attr:`~effmass.analysis.Segment.dk_bohr`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate second order derivatives. Defaults to 100 points evenly distributed across the whole segment.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
        
        Returns:
            array(float). 1d array containing the conductivity effective mass (in units of electron rest mass) evaluated at the points specified in dk.
        """
        dedk, d2edk2 = self.calc_poly_derivatives(polyfit_order=polyfit_order,dk=dk,polyfit_weighting=polyfit_weighting)
        conductivity_mass = [( 1 / x) for x in d2edk2]
        return conductivity_mass

    def calc_transport_effmass(self, polyfit_order=6,dk=[],polyfit_weighting=False):
        r"""
        Calculates the transport mass (:math:`\frac{k}{\delta E \delta k}`), evaluated at :attr:`~effmass.analysis.Segment.dk_bohr` .

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            dk (array, optional): distance (bohr :math:^{-1}) from extrema in reciprocal space at which to evaluate first order derivatives. Defaults to 100 points evenly distributed across the whole segment.
            polyfit_weighting (bool, optional): If True, polyfit will be weighted according to occupation of eigenstates. If False, no weighting will be used.
           
        Returns:
            array(float). 1d array containing the transport effective mass (in units of electron rest mass) evaluated at the points specified in dk.
        """
        dedk, d2edk2= self.calc_poly_derivatives(polyfit_order=polyfit_order,dk=dk,polyfit_weighting=polyfit_weighting)
        # dk=0 for first term gives discontinuity
        if dk == []:
            dk = np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)
        transport_mass = [( x / y) for x, y in zip(dk[1:], dedk[1:])]
        transport_mass = np.append(np.array([np.nan]), transport_mass)
        return transport_mass

    def calc_alpha(self, polyfit_order=6,truncate=True):
        r"""
        The transport mass (:math:`\frac{k}{\delta E \delta k}`) is calculated as a function of :attr:`~effmass.analysis.Segment.dk_bohr` and fitted to a straight line. The gradient of this line determines the alpha parameter which is used in the kane dispersion.
        
        See :meth:`effmass.analysis.Segment.calc_transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            truncate (bool, optional): If True, data only up to :meth:`effmass.analysis.Segment.explosion_index` is used. If False, alpha is calculated using data for the whole segment. Defaults to True.

        Returns:
            float: The alpha parameter (hartree :math:`^{-1}`).
        """
        if truncate:
            self.check_kanefit_points(polyfit_order=polyfit_order)
        transport_mass = self.calc_transport_effmass(polyfit_order=polyfit_order,dk=self.dk_bohr,polyfit_weighting=False)
        if truncate:
            idx = self.explosion_index(polyfit_order=polyfit_order)
            gradient, intercept = np.polyfit(self.dE_hartree[1:idx+1], transport_mass[1:idx+1],1)
        else:
            gradient, intercept = np.polyfit(self.dE_hartree[1:], transport_mass[1:],1)
        alpha = np.divide(gradient,2*intercept)
        return alpha          
                                  
    def calc_kane_mass_bandedge(self, polyfit_order=6,truncate=True):
        r"""
        The transport mass (:math:`\frac{1}{\delta E \delta k}`) is calculated as a function of :attr:`~effmass.analysis.Segment.dk_bohr` and fitted to a straight line. The intercept of this line with the y-axis gives a transport mass at bandedge which is used as a parameter in the kane dispersion.
        
        See :meth:`effmass.analysis.Segment.calc_transport_effmass`.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.
            truncate (bool, optional): If True, data only up to :meth:`effmass.analysis.Segment.explosion_index` is used. If False, alpha is calculated using data for the whole segment. Defaults to True.
        
        Returns:
            float: transport mass at bandedge (in units of electron mass).
        """
        if truncate:
            self.check_kanefit_points(polyfit_order=polyfit_order)
        transport_mass = self.calc_transport_effmass(polyfit_order=polyfit_order,dk=self.dk_bohr,polyfit_weighting=False)
        if truncate:
            idx = self.explosion_index(polyfit_order=polyfit_order)
            gradient, intercept = np.polyfit(self.dE_hartree[1:idx+1], transport_mass[1:idx+1],1)
        else:
            gradient, intercept = np.polyfit(self.dE_hartree[1:], transport_mass[1:],1)
        return intercept   # intercept is the band edge transport mass

    def calc_finite_difference_effmass(self):
        """
        The bandedge curvature is calculated using a second order forward finite difference method. This is then inverted to give an effective mass.

        Args:
            None

        Returns:
            float: Bandedge effective mass from finite difference (in units of electron mass).
        """
        #x =  ((32*self.dE_hartree[1]) - (2*self.dE_hartree[2]))/ (12*self.dk_bohr[1]*self.dk_bohr[1])
        x = (self.dE_hartree[2]-(2*self.dE_hartree[1]))/ (self.dk_bohr[1]**2)
        mass = 1/x
        return mass

    def return_finite_difference_fit(self):
        r"""
        Calculates the curvature at bandedge using a finite difference method and then evaluates the corresponding quadratic dispersion along :attr:`~effmass.analysis.Segment.dk_bohr`.

        See :meth:`effmass.analysis.Segment.calc_finite_difference_effmass`.

        Args:
            None

        Returns:
            list(float): list containing energies (hartree). The energies are calculated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr` using the quadratic approximation. 
        """
        m_bandedge = self.calc_finite_difference_effmass()
        values = [((k**2)/(2*m_bandedge)) for k in np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)]
        return values

    def calc_five_point_effmass(self):
        """
        Calculates a polyfit using only 3 DFT points (symmetry makes 5).

        Args:
            None

        Returns:
            float: Bandedge effective mass from polynomial fit (in units of electron mass).

        Notes:
            no weighting is used.
            function primarily for comparison against a weighted fit across multiple points.
        """

        sym_dk = np.array([-self.dk_bohr[2],-self.dk_bohr[1],self.dk_bohr[0],self.dk_bohr[1],self.dk_bohr[2]])
        sym_dE = np.array([self.dE_hartree[2],self.dE_hartree[1],self.dE_hartree[0],self.dE_hartree[1],self.dE_hartree[2]])
        coeff_matrix = np.vstack([sym_dk ** 2]).T 
        coeff = np.linalg.lstsq(coeff_matrix, sym_dE)[0]
        #print (coeff_linalg)
        #coeff = np.polyfit(sym_dk,sym_dE,2)
        #print (coeff)
        mass = 1/(2*coeff[0])
        return mass

    def return_five_point_fit(self):
        """ 
        Calculates the curvature at band edge using a polynomial fit and then evaluates the corresponding quadratic dispersion along :attr:`~effmass.analysis.Segment.dk_bohr`. 
        
        Args:
            None

        Returns:
            list(float): list containing energies (hartree). The energies are calculated at 100 points evenly distributed across :attr:`~effmass.analysis.Segment.dk_bohr` using the quadratic approximation. 
        """
        m_bandedge = self.calc_five_point_effmass()
        values = [((k**2)/(2*m_bandedge)) for k in np.linspace(self.dk_bohr[0],self.dk_bohr[-1],100)]
        return values

    def explosion_index(self, polyfit_order=6):
        r"""
        This will find the index at which there is a change in sign of the second derivative. 

        In the region of this point the first derivative will pass through zero and so the transport mass (:math:`\frac{1}{\delta E \delta k}`) will explode.

        Args:
            polyfit_order (int, optional): order of polynomial used to approximate the dispersion. Defaults to 6.

        Notes:
            This marks the point at which the Kane dispersion is definitely not valid, although it may be that the Kane dispersion is a poor approximation prior to this.
        """
        dedk, d2edk2 = self.calc_poly_derivatives(polyfit_order=polyfit_order, dk=self.dk_bohr,polyfit_weighting=False)
        sign = np.sign(d2edk2)
        signchange = ((np.roll(sign,1)-sign) !=0).astype(int)
        signchange[0] = 0
        if 1 in signchange:
            cutoff = np.where(signchange==1)[0][0]
        else:
            cutoff = len(self.dE_eV)-1
        return cutoff
        
      


        




            


