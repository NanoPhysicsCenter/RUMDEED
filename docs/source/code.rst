Code description
================

The code

Verlet
------
For particle advancement the code uses the Velocity Verlet algorithm. The code for it can be found in **mod_verlet.F90** and the main subroutine is:

.. f:autosubroutine:: velocity_verlet

The updated positions are calculated with:

.. math::
    x_{n+1} = x_n + v_n\Delta t + \frac{1}{2}a_n \Delta t^2

.. f:autosubroutine:: Update_ElecHole_Position

Then the force on each particle is calculate using Coulomb's law, after which the velocity is
update using:

.. math::
    v_{n+1} = v_n + \frac{a_n+a_{n+1}}{2} \Delta t^2

.. f:autosubroutine:: Update_Velocity

Photo emission
--------------
Photo


Field emission
--------------
Field emission

.. math::
    J = \frac{a}{\phi t^2(l)}F^2 exp(-\nu(l)b\phi^{3/2}/F)

See :cite:p:`Forbes08112007` for nu and t

.. math::
    \nu(l) = 1 - l + \frac{1}{6}l \ln(l)

.. math::
    t(l) = 1 + l\left( \frac{1}{9} - \frac{1}{18}\ln(l) \right)

.. math::
    l = \frac{F}{F_\phi} = \frac{e^3}{4\pi\epsilon_0} \frac{F}{\phi^2}

If \(\phi\) is in eV and \(F\) in V/m then

.. math::
  l = \frac{e}{4\pi\epsilon_0} \frac{F}{\phi^2}


Thermal-field emission
----------------------
Thermo-Field
