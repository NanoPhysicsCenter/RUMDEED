Code description
================

Below is a short description of the theory and equations used in the code.

Verlet
------
For particle advancement the code uses a variant of the Velocity Verlet algorithm :cite:p:`PhysRev.159.98` called Beeman's algorithm :cite:`LEVITT1983617,BEEMAN1976130`. The code for it can be found in **mod_verlet.F90** and the main subroutine is called
``Velocity_Verlet``.

The updated positions are calculated with:

.. math::
    \mathbf{x}_{n}(t + \Delta t) = \mathbf{x}_n(t) + \mathbf{v}_n(t)\Delta t + \frac{1}{6}\left(4\mathbf{a}_n(t) - \mathbf{a}_n(t - \Delta t) \right) \Delta t^2

See the ``Update_ElecHole_Position`` subroutine.

Then the force on each particle is calculate using Coulomb's law,

.. math::
    \mathbf{a}_n(t + \Delta t) = \frac{q_n}{m_n} \left( E_0(\mathbf{x}_n(t+\Delta t))
    + \frac{q_n}{4\pi \epsilon_0} \sum_{k\neq n} q_k\frac{\hat{\mathbf{r}}_{nk}}{\mathbf{r}_{nk}^2}
    + \frac{q_n}{4\pi\epsilon_0}\sum_j Q_j\frac{\hat{\mathbf{r}}_{nj}}{\mathbf{r}_{nj}^2} \right)

The force calculation takes into account the vaccum field :math:`E_0` and the field from other particles. Image charges are also taken into account.
See the ``Calculate_Acceleration_Particles`` subroutine and also the ``Calc_Field_at`` subroutine.
After the force calculate the velocity is update using:

.. math::
    \mathbf{v}_{n}(t+\Delta t) = \mathbf{v}_n(t) + \frac{1}{6}\left( 2\mathbf{a}_{n}(t + \Delta t) + 5\mathbf{a}_n(t) - \mathbf{a}_n(t - \Delta t) \right)\Delta t

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
The field emission code uses the Fowler-Nordheim equation to calculate the total number of electrons to emit.
The code for this can be found in **mod_field_emission_v2.F90**.

The Fowler-Nordheim equation :cite:p:`Forbes08112007` is

.. math::
    J = \frac{a_{FN}}{\phi t^2(\ell)}F^2 exp(-\nu(\ell)b_{FN}\phi^{3/2}/F)

where :math:`a_{FN}` and :math:`b_{FN}` are the first and second Fowler-Nordheim constants respectively.
The functions :math:`\nu(\ell)` and :math:`t(\ell)` arise due to image charge effects and are given by

.. math::
    \nu(\ell) = 1 - \ell + \frac{1}{6}\ell \ln(\ell)

.. math::
    t(\ell) = 1 + \ell\left( \frac{1}{9} - \frac{1}{18}\ln(\ell) \right)

.. math::
    \ell = \frac{F}{F_\phi} = \frac{e^3}{4\pi\epsilon_0} \frac{F}{\phi^2}

If \(\phi\) is in eV and \(F\) in V/m then

.. math::
  \ell = \frac{e}{4\pi\epsilon_0} \frac{F}{\phi^2}

.. _field-tip:

Tip stuff


.. _thermal-field:

Thermal-field emission
----------------------
Thermo-Field

.. index:: Verlet, Beeman, code, Fowler-Nordheim