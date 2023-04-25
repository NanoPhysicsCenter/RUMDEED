Code description
================

The code

Verlet
------
For particle advancement the code uses the Velocity Verlet algorithm. The code for it can be found in **mod_verlet.F90** and the main subroutine is called
``velocity_verlet``.

The updated positions are calculated with:

.. math::
    x_{n+1} = x_n + v_n\Delta t + \frac{1}{2}a_n \Delta t^2

See the ``Update_ElecHole_Position`` subroutine.

Then the force on each particle is calculate using Coulomb's law, after which the velocity is
update using:

.. math::
    v_{n+1} = v_n + \frac{a_n+a_{n+1}}{2} \Delta t^2

See the ``Update_Velocity`` subroutine.

.. _photo:

Photo emission
--------------
Photoemission Branch
    
* Photoemission based on Workfunction pattern - see mod_workfunction as well

    * Checkerboard
    * Random generated pattern - WIP
    * Emission tip (Works but probably not as intended)
    * Circle (Should work, exhaustive tests not done)

* Velocity distribution is Gaussian
    * Velocity has Newtonian calculations

.. math::
    v_z = \sqrt{ \frac{ 2( \hbar \omega - \phi ) * q_{0} }{ m_{0} }}

* Output pulse can be Gaussian
    * Quantum efficiency is controlled by amplitude modulation
    * Pulse repetition is possible - WIP

.. math::
    b = \frac{1}{ 2 \pi \sigma^2}

.. math::
    f(step) = A exp{ - b  ( step - \mu )^2 }

.. math::
    QE \equiv \frac{\hbar \omega}{ q_0 } \frac{ J }{ I_{\omega } } \propto ( \hbar \omega - \Phi)^2

Quantum efficiency and light intensity are both variables in the amplitude so care should be taken when interpreting the input A(mplitude) of the normal distribution.

* Input file for laser and pulse

    * Input file is same as main
    * laser file is photoemission specific

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
