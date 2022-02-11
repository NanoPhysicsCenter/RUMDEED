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

.. _photo:

Photo emission
--------------
Photoemission Branch
    
* Photoemission based on Workfunction pattern - see mod_workfunction as well

    * Checkerboard
    * Random generated pattern - WIP
    * Emission tip
    * Circle - WIP

* Velocity distibution is Gaussian

    * Velocity has Newtonian calculations

* Output pulse can be Gaussian

    * Quantum efficiency is controlled by amplitude modulation
    * Pulse repetition is possible - WIP

* Input file for laser and pulse

    * Input file is same as main
    * laser file is photo_emission specific

Input file

.. code-block:: console

    &
    INPUT
    V_S = 1.0d1,
    BOX_DIM = 0.0d0, 0.0d0, 2500.0d0,
    TIME_STEP = 0.25d-3,
    STEPS = 60000,
    EMISSION_MODE = 1,
    NREMIT = 1,
    IMAGE_CHARGE = .TRUE.,

    EMITTERS_DIM(1:3, 1) = 500.0d0, 500.0d0, 0.0d0,
    EMITTERS_POS(1:3, 1) = -250.0d0, -250.0d0, 0.0d0,
    EMITTERS_TYPE(1) = 2,
    EMITTERS_DELAY(1) = 0,

    /


Laser file

.. code-block:: console

    1 2 2
    4.7 0.02
    10000 1000 5

Input Warning
+++++++++++++

    The header/first line sets parameters;
    The first number enables Gaussian electron emission pulse, 1 = on, 2 = off.
    Second number selects type of laser input, 1 for fixed photon energy, 2 for Poisson distributed photon energy.
    Third number picks velocity profile for electrons, 1 being zero initial velocity, 2 for work function dependant inital velocity.
    
    Second line is laser (photon) energy and variation, first being the laser 'mean' energy level in electronVolts (eV) and second being standard deviation of the laser (in eV's as well). 
    This is normal distribution with Box-Muller method.
    For work function dependant initial velocity the energy is compared to the work function with the excess making way for Newtonian velocity given to the electrons.

    Third line is gauss pulse parameters, center (mu), width (sigma) and A(mplitude) of the pulse. 
    The gaussian pulse is simulated with output restriction of electrons according to normal distribution.
    This should in theory simulate the Quantum Efficiency and Intensity via amplitude modulation.
        

.. _field:

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

.. _field-tip:

Tip stuff


.. _thermal-field:

Thermal-field emission
----------------------
Thermo-Field
