==============
Implementation
==============

Methods for calculating the curvature effective mass
====================================================

Three methods are used to calculate the curvature effective mass (Eqn.\ 1 in the main text).

- **Finite difference**

  We use a three point forward finite difference equation to calculate the curvature at point i:

  .. math::

    \frac{\partial^2E}{\partial k^2} = \frac{E_{i+2} - 2E_{i+1} + E_{i}}{\left|k_{i+1} - k_i\right|},

  where :math:`E_i` is the energy eigenvalue at position :math:`k_i` in reciprocal space. 


- **Unweighted least-squares fit**

  To obtain estimates for the coefficient of a parabolic dispersion

  .. math::

    E = ck^2,


  we use the least-squares method as implemented in the NumPy Python library to minimise the summed square of residuals

  .. math::

    \sum^{5}_{i=1}(ck_i^2 - E_i)^2.


  We fit to five points; three points from the DFT-calculated dispersion plus two from the symmetry of the dispersion (:math:`E(k)=E(-k)`).
  

- **Weighted least-squares fit**

  To obtain estimates for the coefficients of the dispersion

  .. math::

    E = {c}k^2,


  we use the least-squares method as implemented in the NumPy Python library to minimise the summed square of residuals

  .. math::

    \sum^{n}_{i=1}W_i(ck_i^2 - E_i)^2.


  The summation is over all points up to an energy of :math:`0.25\,eV`, including points generated from the symmetry of the dispersion, :math:`E(k)=E(-k)`.
  :math:`W_i` is given by

  .. math::

    W_i(E_i,T) = \frac{1}{\exp\left(\frac{E_i-E_f}{k_BT}\right)+1}.


