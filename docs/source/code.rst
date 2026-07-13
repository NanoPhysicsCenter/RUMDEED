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

Then the force on each particle is calculated using Coulomb's law,

.. math::
    \mathbf{a}_n(t + \Delta t) = \frac{q_n}{m_n} \left( E_0(\mathbf{x}_n(t+\Delta t))
    + \frac{q_n}{4\pi \epsilon_0} \sum_{k\neq n} q_k\frac{\hat{\mathbf{r}}_{nk}}{\mathbf{r}_{nk}^2}
    + \frac{q_n}{4\pi\epsilon_0}\sum_j Q_j\frac{\hat{\mathbf{r}}_{nj}}{\mathbf{r}_{nj}^2} \right)

The force calculation takes into account the vacuum field :math:`E_0` and the field from other particles. Image charges are also taken into account.
See the ``Calculate_Acceleration_Particles`` subroutine and also the ``Calc_Field_at`` subroutine.
After the force calculation the velocity is updated using:

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
    v_z = \sqrt{ \frac{ 2( \hbar \omega - \phi ) q_{0} }{ m_{0} }}

* Output pulse can be Gaussian
    * Quantum efficiency is controlled by amplitude modulation
    * Pulse repetition is possible - WIP

.. math::
    b = \frac{1}{ 2 \pi \sigma^2}

.. math::
    f(step) = A \exp\left( - b ( step - \mu )^2 \right)

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
    J = \frac{a_{FN}}{\phi t^2(\ell)}F^2 \exp(-\nu(\ell)b_{FN}\phi^{3/2}/F)

where :math:`a_{FN}` and :math:`b_{FN}` are the first and second Fowler-Nordheim constants respectively.
The functions :math:`\nu(\ell)` and :math:`t(\ell)` arise due to image charge effects and are given by

.. math::
    \nu(\ell) = 1 - \ell + \frac{1}{6}\ell \ln(\ell)

.. math::
    t(\ell) = 1 + \ell\left( \frac{1}{9} - \frac{1}{18}\ln(\ell) \right)

.. math::
    \ell = \frac{F}{F_\phi} = \frac{e^3}{4\pi\epsilon_0} \frac{F}{\phi^2}

If :math:`\phi` is in eV and :math:`F` in V/m then

.. math::
  \ell = \frac{e}{4\pi\epsilon_0} \frac{F}{\phi^2}

.. _field-tip:

Field emission from a hyperboloid tip
-------------------------------------
The code for field emission from a hyperboloid tip can be found in **mod_emission_tip.f90**.
The work function of the tip is set with the ``w_theta`` variable in that file.

The tip is modeled as a hyperboloid of revolution, which is described using prolate spheroidal coordinates.
With :math:`\xi = \cosh{\mu}` and :math:`\eta = \cos{\nu}` the coordinates are defined as

.. math::
    \begin{split}
    x &= a \sqrt{\xi^2-1}\sqrt{1-\eta^2} \cos{\phi}\, ,\\
    y &= a \sqrt{\xi^2-1}\sqrt{1-\eta^2} \sin{\phi}\, ,\\
    z &= a \xi \eta\, ,
    \end{split}

where :math:`a` is the focal distance of the hyperboloid. The reverse relations are

.. math::
    \begin{split}
    \xi &= \frac{1}{2a} \left( \sqrt{x^2 + y^2 + (z+a)^2} + \sqrt{x^2 + y^2 + (z-a)^2} \right)\, ,\\
    \eta &= \frac{1}{2a} \left( \sqrt{x^2 + y^2 + (z+a)^2} - \sqrt{x^2 + y^2 + (z-a)^2} \right)\, ,\\
    \phi &= \arctan{\frac{y}{x}}\, .
    \end{split}

The surface of the tip is the hyperboloid :math:`\eta = \eta_1`, which is held at :math:`V = 0`,
and the anode is the plane :math:`\eta_2 = 0` at :math:`V = V_0`.
The potential between them is :cite:p:`pan:2151`

.. math::
    V(\eta) = V_0 \frac{\ln{\left[ \frac{1 + \eta_1}{1-\eta_1}\frac{1-\eta}{1+\eta} \right]}}{\ln{\left[ \frac{1+\eta_1}{1-\eta_1}\frac{1-\eta_2}{1+\eta_2} \right]}}\, .

Note that the boundary conditions have been swapped compared to the reference.
The vacuum electric field is then

.. math::
    \vec{E} = - \nabla V(\eta) = \frac{2V_0}{a} \frac{1}{\xi^2-\eta^2} \frac{1}{\ln \left[ \frac{1+\eta_1}{1-\eta_1} \frac{1-\eta_2}{1+\eta_2} \right]}
    \begin{pmatrix}
    -\eta \sqrt{\frac{\xi^2-1}{1-\eta^2}} \cos{\phi}\\
    -\eta \sqrt{\frac{\xi^2-1}{1-\eta^2}} \sin{\phi}\\
    \xi
    \end{pmatrix}\, ,

with the magnitude

.. math::
    |\vec{E}| = \frac{2V_0}{a} \frac{1}{\sqrt{\xi^2-\eta^2}\sqrt{1-\eta^2}} \frac{1}{\ln \left[ \frac{1+\eta_1}{1-\eta_1} \frac{1-\eta_2}{1+\eta_2} \right]}\, .

At the apex of the tip, where :math:`\eta = \eta_1` and :math:`\xi = 1`, the field points in the z direction,

.. math::
    E_z = \frac{2V_0}{a} \frac{1}{1 - \eta_1^2} \frac{1}{\ln \left[ \frac{1+\eta_1}{1-\eta_1} \right]}\, .

The shape of the tip is set in the input file with the parameter :ref:`EMITTERS_DIM <run>`,
which gives the distance :math:`d` from the apex of the tip to the anode, the base radius :math:`R` and the height :math:`h` of the tip.
These determine the coordinate parameters through

.. math::
    R = a \sqrt{\xi_{max}^2-1}\sqrt{1-\eta_1^2}\, , \quad
    \eta_1 = -\frac{d}{a}\, , \quad
    \xi_{max} = \frac{h}{d} + 1\, , \quad
    a = \sqrt{\frac{d^2R^2}{h^2 + 2dh} + d^2}\, ,

which keep the shape of the tip constant for all values of :math:`d`.
The radius of curvature of the surface is

.. math::
    R_c = \left| \frac{a}{\eta_1} \frac{(\xi^2 - \eta_1^2)^{\frac{3}{2}}}{\sqrt{1-\eta_1^2}} \right| \, ,

which at the apex, where :math:`\xi = 1` and :math:`\eta_1 = -\frac{d}{a}`, reduces to

.. math::
    R_c = \frac{a^2}{d} - d\, .

The surface area of the tip between :math:`\xi_1` and :math:`\xi_2`, used when integrating the emission current over the surface, is

.. math::
    A = \frac{a^2}{2} \sqrt{1-\eta_1^2}\, (\phi_2 - \phi_1) \left[ \xi \sqrt{\xi^2 - \eta_1^2} - \eta_1^2 \ln\left(\xi + \sqrt{\xi^2 - \eta_1^2}\right) \right]_{\xi_1}^{\xi_2}\, .

For the image charge interaction the tip is approximated by a sphere.
A charged particle at a distance :math:`y` from the center of a sphere with radius :math:`r` has an image charge partner at a distance

.. math::
    y^\prime = \frac{r^2}{y}

from the center of the sphere, with the charge

.. math::
    q^\prime = -\frac{r}{y}q\, ,

where :math:`q` is the charge of the particle at :math:`y`.

.. _thermal-field:

Thermal-field emission
----------------------
The code for planar thermal-field emission can be found in **mod_field_thermo_emission.F90** :cite:p:`PhysRevApplied.15.014040`.

.. index:: Verlet, Beeman, code, Fowler-Nordheim