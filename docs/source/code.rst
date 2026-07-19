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

.. _twotimestep:

Two-time-step integration (collision runs)
------------------------------------------
In collision runs (see :ref:`COLLISION_MODE <run>`) the system contains ions and a neutral background
gas in addition to the electrons. The ions are more than four orders of magnitude heavier than the
electrons and barely move during an electron time step, so integrating them with the electron
:math:`\Delta t` wastes work. Setting :ref:`ATOM_TIME_INTERVAL <run>` to an integer :math:`n > 0`
activates a two-time-step (multiple-time-step) scheme in **mod_verlet.F90**:

* **Electrons** take a Beeman step of :math:`\Delta t` on every time step, exactly as in the
  standard integrator. They feel the ions at their frozen positions.
* **Ions** take a Beeman step of :math:`n \Delta t` on every :math:`n`-th time step
  (``Update_Species_Position`` / ``Update_Species_Velocity`` with the step length scaled by
  :math:`n`). Between those steps their positions are frozen.
* **Atoms** are never advanced; they act as a static background for the collision sampling.

The forces on an ion are split by time scale, following the impulse multiple-time-step method
(the same idea as r-RESPA :cite:p:`doi:10.1063/1.463137`). The *slow* forces --- ion--ion pairs,
their image charges and the vacuum field --- are evaluated once per ion step
(``Update_Ion_Acceleration``) and integrated with the Beeman formulas above with
:math:`\Delta t \rightarrow n\Delta t`. The *fast* electron--ion pair forces change on the electron
time scale, so they are not accumulated into the Beeman acceleration history. Instead each electron
step applies them to the ion velocity directly as an impulse,

.. math::
    \mathbf{v}_\text{ion} \leftarrow \mathbf{v}_\text{ion} + \mathbf{a}_\text{fast}(t_k)\,\Delta t ,

(see the fold at the end of ``Update_Elec_Acceleration``). Over one ion window the impulses sum to
:math:`\sum_k \mathbf{a}_\text{fast}(t_k)\Delta t \approx \int \mathbf{a}_\text{fast}\,dt`, so the
momentum transferred from the electron cloud to the ion is captured at full time resolution while
the expensive ion--ion loop only runs every :math:`n`-th step. Accumulating the fast samples into
the Beeman acceleration instead would weight them :math:`n` times too strongly, since the ion update
treats that array as a single acceleration and scales it by :math:`n\Delta t` (this was a bug in
older versions, fixed July 2026; see below).

Details worth knowing:

* **Newly created ions.** An ion created by ionization between two ion steps has only existed for
  part of the window, so its first update caps the step length at
  :math:`\min(n,\, \text{step} - \text{birth step} + 1)` steps.
* **Ramo current.** The induced (Shockley--Ramo) current :math:`q\,\mathbf{v}\cdot\hat{\mathbf{E}}`
  of an ion is recorded on *every* time step --- on the ion steps by the Beeman velocity update and
  in between by ``Update_Ion_Ramo_Current`` --- so the current written to **ramo_current.dt**
  integrates to the transported ion charge. The same holds for the ion velocity averages in that
  file. Output produced before the July 2026 fixes under-reports the ion current and the
  between-step ion velocity averages by a factor of :math:`n`.
* **Equivalence at** :math:`n = 1`. With ``ATOM_TIME_INTERVAL = 1`` every step is an ion step and
  the single fast-force sample is exactly the Beeman acceleration, so the two-time-step integrator
  reproduces the standard one bitwise. This contract is asserted by ``Test_TTS_Equivalence`` in
  **mod_tests.F90**; the staggered cadence (:math:`n > 1`), the impulse treatment and the per-step
  Ramo recording are asserted by ``Test_TTS_Staggered``.
* **Accuracy.** The impulse channel carries an :math:`O(\Delta t)` half-step phase offset relative
  to the Beeman weighting of the reference integrator (about :math:`0.5/N_\text{steps}` relative on
  velocity after :math:`N_\text{steps}` steps), and the fast impulses gathered during a window only
  affect the ion *drift* from the next window on, an :math:`O((n\Delta t)^2)` per-window effect.
  Choose :math:`n` so that the ions move a negligible fraction of the inter-particle spacing per
  window; the electron--ion forces themselves are sampled at full resolution regardless of
  :math:`n`.

.. _collisions:

Electron collisions with the background gas
-------------------------------------------
In collision runs the system contains a neutral background gas of nitrogen (N\ :sub:`2`) in addition
to the emitted electrons. Energetic electrons ionize the gas, creating an electron–ion pair, and the
resulting ions can later recapture an electron by recombination. The code lives in
**mod_collisions.F90**; the top-level routine ``Do_Electron_Atom_Collisions`` is called once per time
step from **mod_verlet.F90** and dispatches on :ref:`COLLISION_MODE <run>`:

.. list-table::
   :header-rows: 1
   :widths: 10 45 45

   * - Mode
     - Ionization
     - Recombination
   * - 1
     - Continuous
     - —
   * - 2
     - Continuous
     - Discrete
   * - 3
     - Discrete
     - —
   * - 4
     - Discrete
     - Discrete

The *continuous* ionization model is due to Kristinn Torfason. The *discrete* ionization model and the
recombination model are due to Arnar Jónsson.

Background gas and cross sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The gas is treated as stationary and uniform, with number density set from the ideal gas law from the
input temperature :ref:`T_TEMP <run>` and pressure :ref:`P_ABS <run>`,

.. math::
    n_d = \frac{P}{k_B T}\, .

Two electron-impact cross sections for N\ :sub:`2` enter the models: the total (momentum-transfer)
cross section :math:`\sigma_\text{tot}(E)` and the ionization cross section :math:`\sigma_\text{ion}(E)`,
both a function of the electron kinetic energy :math:`E` (in eV). They are tabulated in the files
``N2-tot-cross.txt`` and ``N2-ion-cross.txt`` (read by ``Read_Cross_Section``). Below the fitted range
the value is linearly interpolated from the table (clamped to the tabulated energy range); above it a
double-exponential fit :math:`\sigma(E) = (a\, e^{bE} + c\, e^{dE})\times 10^{-20}\ \text{m}^2` is used
(``Find_Cross_tot_data``, ``Find_Cross_ion_data``). For cross-section lookups the energy is capped at
5000 eV.

Each electron caches its energy and the associated cross-section quantities in
``Update_Collision_Data``, refreshed for all electrons at the start of every collision step. This
routine also converts each cross section to an equivalent geometric interaction radius by treating it
as the area of a disk, :math:`\sigma = \pi r^2`, so

.. math::
    r_\text{ion} = \sqrt{\sigma_\text{ion}/\pi}\, , \qquad
    r_\text{rec} = \sqrt{\sigma_\text{rec}/\pi}\, ,

which the discrete models below use as collision radii.

Continuous ionization
~~~~~~~~~~~~~~~~~~~~~~~
*Model by Kristinn Torfason* (``Do_Continuous_Ionization_ots`` / ``_tts``).

In the continuous model the gas is not represented by explicit particles; each electron instead
samples a collision probability from the local mean free path. An electron that moved a distance
:math:`d = |\mathbf{x}(t) - \mathbf{x}(t-\Delta t)|` in the last step, and whose energy exceeds the
binding energy :math:`E > E_\text{bind}`, has mean free path and collision probability

.. math::
    \lambda = \frac{1}{n_d\, \sigma_\text{tot}(E)}\, , \qquad
    P_\text{coll} = \frac{d}{\lambda}\, ,

valid while :math:`d \ll \lambda`. A uniform random number decides whether a collision occurs. If it
does, the collision is ionizing with the conditional probability

.. math::
    P_\text{ion} = \frac{\sigma_\text{ion}(E)}{\sigma_\text{tot}(E)}\, ,

otherwise the collision is treated as elastic and leaves the electron unchanged (only the collision is
counted).

An ionizing collision converts an incident electron of energy :math:`E_1` into two outgoing electrons.
The binding energy is paid first and the remainder is split randomly between the incident (scattered)
electron and the newly ejected electron,

.. math::
    E_2 = E_1 - E_\text{bind}\, , \qquad
    E_\text{scat} = \zeta E_2\, , \qquad
    E_\text{ejec} = (1-\zeta) E_2\, , \qquad \zeta \sim U(0,1)\, ,

so that energy is conserved. Both electrons are given new directions (see
:ref:`collision-angles`) and speeds :math:`v = \sqrt{2 q_0 E/m_0}`. The ejected electron is placed at a
small random offset from the collision point, and a singly charged ion of mass
:math:`m_{\text{N}_2^+}` is created at rest nearby with a finite lifetime :ref:`ION_LIFE_TIME <run>`.

Discrete ionization
~~~~~~~~~~~~~~~~~~~~~
*Model by Arnar Jónsson* (``Do_Discrete_Ionization_ots`` / ``_tts``).

In the discrete model the gas is represented by explicit neutral-atom particles, placed uniformly in a
cylinder over the emitter at the start of the run by ``Place_Atoms_Cylinder`` (a fraction
:ref:`ION_ATOM_RATIO <run>` of them start already ionized). Ionization is then a geometric event: an
electron ionizes an atom if its trajectory over the time step passes within the ionization radius
:math:`r_\text{ion}` of the atom, provided the electron has enough energy, :math:`E \ge E_\text{bind}`.

The electron is advanced over the step with the Beeman trajectory
:math:`\mathbf{x}(t) = \mathbf{x}_0 + \mathbf{v} t + \tfrac{1}{2}\mathbf{a} t^2`. Writing
:math:`\mathbf{r}_0 = \mathbf{x}_0 - \mathbf{x}_\text{atom}`, the crossing condition
:math:`|\mathbf{x}(t) - \mathbf{x}_\text{atom}|^2 = r_\text{ion}^2` is the quartic

.. math::
    \tfrac{1}{4}|\mathbf{a}|^2\, t^4
    + (\mathbf{v}\cdot\mathbf{a})\, t^3
    + \left(|\mathbf{v}|^2 + \mathbf{r}_0\cdot\mathbf{a}\right) t^2
    + 2(\mathbf{r}_0\cdot\mathbf{v})\, t
    + \left(|\mathbf{r}_0|^2 - r_\text{ion}^2\right) = 0\, ,

solved for its roots by ``SolvePolynomial`` in **mod_polynomialroots.F90**. If the electron already
starts inside the radius the event fires immediately; otherwise the earliest real root
:math:`t \in (0, \Delta t]` is the crossing time. The energy split and the angular distributions of the
two outgoing electrons are the same as in the continuous model. On ionization the atom is removed and
replaced by an ion at the atom's position, and the ejected electron is created nearby, so that charge
and particle number are accounted for exactly.

Recombination
~~~~~~~~~~~~~
*Model by Arnar Jónsson* (``Do_Discrete_Recombination_ots`` / ``_tts``).

Recombination is the reverse process: a free electron is captured by an ion. It uses the same geometric
crossing test as discrete ionization, but with a capture radius :math:`r_\text{rec}` derived from the
Kramers radiative-recombination cross section :cite:p:`kramers1923`
(``Calculate_Kramers_Cross_Section``),

.. math::
    \sigma_\text{rec}(E) = 2.105\times 10^{-26}\,
    \frac{\text{Ry}^2\, Z_\text{eff}^4}
         {N_n E\left(N_n^2 E + \text{Ry}\, Z_\text{eff}^2\right)}\, ,

where :math:`\text{Ry}` is the Rydberg energy, :math:`N_n = 2` the principal quantum number of nitrogen
and :math:`Z_\text{eff}` an effective charge fixed from the binding energy. Because this cross section
diverges as :math:`1/E` at low energy, the energy used to evaluate it is floored at the thermal energy
:math:`k_B T` of the gas, which caps the capture radius for slow electrons.

For every ion–electron pair the electron's trajectory over the step is tested against :math:`r_\text{rec}`
with the same quartic as above. When an electron crosses within the capture radius, the ion and the
electron are both removed. In mode 4 (where the gas is made of explicit atoms) a neutral atom is put
back at the ion's position, so the discrete gas is conserved over an ionize–recombine cycle; in mode 2
(continuous ionization) no atom is returned. Ions are also removed when they reach the end of their
lifetime :ref:`ION_LIFE_TIME <run>`.

.. _collision-angles:

Scattering angles
~~~~~~~~~~~~~~~~~~
Both outgoing electrons of an ionizing collision are given a direction relative to the incident
electron's velocity, drawn by acceptance–rejection from a folded normal distribution
:cite:p:`leone1961folded` of the scattering angle :math:`\theta \in [0^\circ, 180^\circ]`,

.. math::
    f(\theta) = \sqrt{\frac{2}{\pi \sigma^2}}\,
    \exp\!\left(-\frac{\mu^2 + \theta^2}{2\sigma^2}\right) \cosh\!\left(\frac{\mu\theta}{\sigma^2}\right)\, ,

with the envelope constant taken as the true maximum of :math:`f` over the angular domain
(``folded_normal_max``) so the acceptance ratio never exceeds one.

* The **scattered** (incident) electron uses a distribution peaked near the forward direction,
  :math:`\mu = 5^\circ`, :math:`\sigma = 25^\circ` (``Get_Injected_Vec``).
* The **ejected** electron uses an energy-dependent peak angle
  :math:`\theta_\text{max}(T) = a\, T^{b} + c` with :math:`a = -430.5`, :math:`b = -0.5445`,
  :math:`c = 89.32` (clamped for :math:`T < 100` eV) and :math:`\sigma = 48^\circ`
  (``Get_Ejected_Vec``). This form follows the doubly differential cross-section measurements for
  electrons scattered by nitrogen of Nogueira, Eschiapati Ferreira and Barbieri.

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

.. _metropolis-hastings:

Metropolis-Hastings sampling of the emission positions
------------------------------------------------------
In every time step the emission happens in two stages :cite:p:`doi:10.1063/1.4914855`.
First the number of electrons to try to emit, :math:`N`, is found by integrating over the surface of the emitter.
Then each of the :math:`N` candidate electrons is given a position on the surface, drawn with the Metropolis-Hastings algorithm.

For field emission the current density is split into the product

.. math::
    J = S \cdot D\, , \quad
    S = \frac{a_{FN}}{\phi\, t^2(\ell)}F^2\, , \quad
    D = \exp(-\nu(\ell)\, b_{FN}\, \phi^{3/2}/F)\, ,

where :math:`S` is the electron supply and :math:`D` the escape probability (see :ref:`field`).
For the emitted positions to follow the current density, the number :math:`N`, the distribution the chains sample from, and the emission test applied to each candidate must fit together.
The code uses one of two consistent pairings:

* **Supply pairing**: :math:`N = \frac{\Delta t}{e}\int S\, \mathrm{d}A`, the chains sample positions with density proportional to :math:`S`, and each candidate is then emitted with the probability :math:`D` evaluated at its position.
  The emitted positions follow :math:`S \cdot D = J` and the expected number of emitted electrons is :math:`\frac{\Delta t}{e}\int J\, \mathrm{d}A`.
  Used by the planar field emission and by the field emission from the tip.

* **Current-density pairing**: :math:`N` is drawn from a Poisson distribution with mean :math:`\frac{\Delta t}{e}\int J\, \mathrm{d}A`, the chains sample positions proportional to :math:`J`, and every candidate is emitted.
  Used by the planar thermal-field emission, the thermal-field (GTF) emission from the tip, and the field emission from the torus.

Note that the chains must never sample proportional to :math:`D` (or :math:`J`) when the caller also applies the emission test :math:`D`: the escape probability would then be counted twice and the emitted positions would follow :math:`D^2` instead of :math:`J`.
The samplers therefore always target the expression whose integral sets :math:`N`, divided by whatever probability the caller applies afterwards.

**The chains.**
Each candidate electron gets its own chain.
The chain starts at a random position, uniform in the surface parameters, redrawn until the field there is favorable for emission.
A jump proposes a new position by adding a Gaussian step to the surface parameters.
Steps that leave the emitter are reflected at its edges, and angles wrap around, both of which keep the proposal symmetric so the plain Metropolis acceptance rule applies.
The acceptance is evaluated in log space because the current density spans many orders of magnitude,

.. math::
    \alpha = \ln \pi_\mathrm{new} - \ln \pi_\mathrm{cur}\, , \quad
    \text{accept if}\ \ln r \le \alpha\, , \quad r \sim U(0,1)\, ,

where :math:`\pi` is the target density.
A proposal where the field is not favorable for emission has zero target and is always rejected.
When the chain moves in surface parameters :math:`(u, v)` instead of directly on the surface, the target includes the area element of the parametrization, :math:`\pi \propto J\, h_u h_v` with :math:`\mathrm{d}A = h_u h_v\, \mathrm{d}u\, \mathrm{d}v`, so that the stationary distribution is proportional to the current density per unit area.
Constant factors cancel in the Metropolis ratio and are left out of the targets in the code.

**Self-tuning.**
The first quarter of the jumps of a chain (the warmup) use a large fixed step, 10% of the parameter ranges, so that the chain finds the high-current region from any starting point.
The remaining jumps use the adapted step size :math:`\sigma`, which is shared between the chains, held fixed within a chain, and updated once at the end of each chain from that chain's acceptance rate :math:`a`,

.. math::
    \sigma \leftarrow \sigma\, e^{\gamma (a - a^*)}\, ,

clamped between a lower and an upper limit.
A high acceptance rate means the steps are too small, so the step grows, and the other way around.
The target rate is :math:`a^* = 35\%` in all modules, the optimum for a random walk sampler in two dimensions.
The choice is not critical: the step size only affects how fast the chain decorrelates, never the distribution it samples, and the sampling efficiency is close to optimal for acceptance rates anywhere between roughly 20% and 50%.

Planar emitters
~~~~~~~~~~~~~~~
Found in **mod_field_emission_v2.F90**.
Field emission uses the supply pairing: ``Do_Field_Emission_Planar_rectangle`` gets :math:`N` from the surface integral of :math:`S` and emits each candidate with the probability :math:`D` (``Escape_Prob_log``) at the sampled position.
The chains (``Metropolis_Hastings_rectangle_J``, 200 jumps) move directly in :math:`(x, y)` on the emitter, where the area element is trivial, and target

.. math::
    \ln S = 2\ln|F| - 2\ln t(\ell) - \ln \phi

(``Elec_Supply_log``).
A batched version, ``Metropolis_Hastings_rectangle_J_batch``, advances all chains of the time step in lockstep, so the surface field of every chain is evaluated in one batch per jump — a single GPU kernel per jump iteration in OpenACC builds.
It is enabled with the ``MH_BATCH`` flag in the input file and its step-size update is applied once per jump iteration, from the acceptance rate of that iteration across all chains.
The ring emitter (``Metropolis_Hastings_ring_J``) moves in the ring parameters, angle and radius: the angle wraps around the ring, proposals outside the radial limits are rejected, and the target carries the :math:`\ln r` Jacobian of the parametrization.

Planar thermal-field emission (**mod_field_thermo_emission.F90**) uses the current-density pairing: :math:`N` is Poisson distributed with the surface integral of the general thermal-field (GTF) current density as its mean (``Get_Kevin_Jgtf_v2``), the chains (25 jumps, no warmup) target :math:`\ln J_{GTF}`, and every candidate is emitted.

Field emission from 2D materials (**mod_field_emission_2D.f90**, the 2DEG and Dirac models) is a special case of the supply pairing where no chain is needed at all: the electron supply of these models is a material constant times the emitter area, so the supply density over the surface is exactly uniform and the consistent candidate distribution is the uniform one (a Metropolis chain with a constant target accepts every jump and is just a wasted random walk).
``Uniform_Pos_rectangle`` draws each candidate uniformly over the emitter and the caller emits it with the escape probability :math:`D(F)` at that spot, so the emitted positions follow the current density.
A candidate that lands where the field is not favorable is discarded rather than redrawn, because :math:`N` counts the supply of the whole area, so the supply of a blocked spot must be lost instead of redistributed.

Hyperboloid tip
~~~~~~~~~~~~~~~
Found in **mod_emission_tip.f90**.
The chains move in the :math:`(\xi, \phi)` parametrization of the tip surface (see :ref:`field-tip`), where the area element is

.. math::
    h_\xi h_\phi = a^2 \sqrt{(\xi^2 - \eta_1^2)(1 - \eta_1^2)}\, ,

so the targets carry an extra :math:`\frac{1}{2}\ln(\xi^2 - \eta_1^2)` term.
Steps in :math:`\xi` are reflected at the apex :math:`\xi = 1` and at :math:`\xi_{max}`, and :math:`\phi` wraps around the tip.
The chains are 80 jumps long.

* Field emission (``EMITTERS_TYPE = 1``) uses the supply pairing: ``Do_Field_Emission_Tip_OLDCODE`` integrates :math:`S` over the surface with a midpoint rule, the chains (``Metro_algo_tip_v3``) target :math:`\ln(S\, h_\xi h_\phi)` (``Tip_fe_target_log``), and each candidate is emitted with the escape probability at the sampled position.

* Thermal-field (GTF) emission (``EMITTERS_TYPE = 4``) uses the current-density pairing: ``Do_GTF_Emission_Tip`` integrates :math:`J_{GTF}\, h_\xi h_\phi` with Cuba, the chains (``Metro_algo_tip_gtf_v3``) target the same expression (``Tip_gtf_target_log``), and every candidate is emitted.

The older samplers ``Metro_algo_tip_v2``, ``Metro_algo_tip_gtf`` and ``Metro_algo_tip_gtf_v2`` are kept in the file for reference but are superseded: they target the escape probability instead of the quantity that is integrated for :math:`N`, and their step sizes were effectively frozen (see the comments in the code).

Torus
~~~~~
Found in **mod_torus.F90**.
The torus (a looped carbon nanotube) is parametrized by the angle :math:`\phi` along the loop, with radii :math:`R_y` and :math:`R_z`, and the poloidal angle :math:`\theta` around the tube cross-section of radius :math:`\rho`,

.. math::
    \begin{split}
    x &= \rho \sin\theta\, ,\\
    y &= (R_y + \rho\cos\theta) \cos\phi\, ,\\
    z &= (R_z + \rho\cos\theta) \sin\phi\, .
    \end{split}

The area element is :math:`h_\phi h_\theta` with :math:`h_\theta = \rho` and

.. math::
    h_\phi = \sqrt{ (R_y + \rho\cos\theta)^2 \sin^2\phi + (R_z + \rho\cos\theta)^2 \cos^2\phi }\, .

The vacuum field is interpolated from a finite-element mesh with a kd-tree.
Emission uses the current-density pairing: :math:`N` is Poisson distributed with the Cuba integral of :math:`J\, h_\phi h_\theta` over the loop as its mean, and every candidate is emitted.
The chains (``Metropolis_Hastings_Torus_v2``, 500 jumps) are restricted to a :math:`\pm 10^\circ` patch around the top of the loop (:math:`\phi = \pi/2`, :math:`\theta = 0`), where the field is highest and outside of which the current density is negligible, and target :math:`\ln(J\, h_\phi)` (``Torus_target_log``, :math:`h_\theta` is constant).
The original sampler ``Metropolis_Hastings_Torus``, which targets :math:`|F|` without the area element, is kept for reference.

Cylindrical tip
~~~~~~~~~~~~~~~
Found in **mod_cylindrical_tip.F90**.
The emitter is a cylinder of radius :math:`R` and height :math:`h` whose top edge is rounded with a quarter-torus corner of radius :math:`r_c`.
Its surface consists of three sections — the top disk, the corner arc and the side — but the chain does not move in the local coordinates of the sections.
Instead it uses one seamless chart of the whole surface of revolution: the meridian arc length :math:`s`, running from the center of the top (:math:`s = 0`) over the top disk, around the corner arc and down the side, together with the azimuth :math:`\phi`.
A Gaussian step in :math:`(s, \phi)` is then symmetric everywhere, also across the section boundaries, so no jumps between coordinate patches are needed.
The area element of this chart is :math:`\mathrm{d}A = r(s)\, \mathrm{d}s\, \mathrm{d}\phi`, where :math:`r(s)` is the distance from the axis (``Cyl_meridian_radius``),

.. math::
    r(s) =
    \begin{cases}
    s & \text{top disk,} \\
    (R - r_c) + r_c \sin\left( \frac{s - s_1}{r_c} \right) & \text{corner arc,} \\
    R & \text{side,}
    \end{cases}

with :math:`s_1 = R - r_c` the edge of the top disk.
The vacuum field is interpolated from a finite-element mesh with a kd-tree.
Emission uses the current-density pairing: :math:`N` is Poisson distributed with the Cuba integral of :math:`J` over the three sections as its mean, and every candidate is emitted.
The chains (``Metropolis_Hastings_cyl_tip_v2``, 55 jumps) target :math:`\ln(J\, r(s))` (``Cyl_target_log``); :math:`s` is reflected at the center of the top and at the bottom of the side, and :math:`\phi` wraps around the axis.
Because the current density concentrates almost entirely on the corner arc, which is a small fraction of the meridian, the chain starts in a section drawn with the probabilities that the surface integration computes from the per-section current integrals each time step, then uniform along that section's meridian — a starting point only, not a bias, since the chain converges to its target from any start.
The original sampler ``Metropolis_Hastings_cyl_tip`` is kept for reference: it proposes steps in the local coordinates of each section and projects points that cross a boundary onto the neighboring section; those transition maps are not symmetric, the physical step sizes differ by large factors between the sections, and it targets :math:`|F|` without the area elements.

.. _tip-ic:

Image charge on the prolate spheroidal tip
------------------------------------------
The main code approximates the tip by a sphere when it computes the image charge interaction, as described at the end of :ref:`field-tip`.
**Util/Tip-IC** is a standalone utility that instead computes the image charge on the prolate spheroidal tip itself, using the method of :cite:p:`peridier:4888`, and can be used to check how good the spherical approximation is.
It is not part of the simulation; it is built separately with its own makefile in ``Util/Tip-IC/build`` and writes the induced surface charge density on a grid to ``charge_graph.dt``.

The geometry is the same as in :ref:`field-tip`: the tip is the surface :math:`\eta = \eta_1` in prolate spheroidal coordinates, the anode is the plane :math:`\eta_2 = 0`, and a point charge sits at :math:`(\xi_c, \eta_c, \phi_c)`.
The induced charge density is obtained from an integral over the separation constant :math:`\tau` of a sum over the azimuthal index :math:`m`,

.. math::
    \sigma(\xi, \phi) = \int_0^\infty \sum_{m=0}^{30} H_m(\tau, \xi, \phi) \, \mathrm{d}\tau\, ,

Writing :math:`K^m_\tau(u) = P^m_{-\frac{1}{2}+i\tau}(u)` for the conical function, each term separates into a part that depends only on :math:`(m, \tau)` and the fixed source and tip geometry, and a part that depends on the field point,

.. math::
    H_m(\tau, \xi, \phi) = A_m(\tau)\, \cos(m\phi)\, K^m_\tau(\xi)\, ,

with

.. math::
    A_m(\tau) = (-1)^m \left(2 - \delta_{0m}\right) \tau \tanh(\pi\tau)\,
    \frac{\Gamma\!\left(\frac{1}{2} - m + i\tau\right)}{\Gamma\!\left(\frac{1}{2} + m + i\tau\right)}\,
    K^m_\tau(\xi_c)\,
    \frac{K^m_\tau(-\eta_c) - K^m_\tau(\eta_c)}{K^m_\tau(\eta_1) - K^m_\tau(-\eta_1)}\, ,

where the last factor is the boundary term with :math:`\eta_2 = 0` factored out.
Because :math:`A_m` does not depend on the field point and the :math:`\tau` quadrature points form a fixed grid, it is tabulated once by ``Precompute_Hm_Prefactors`` instead of being recomputed for every point of the output grid, with a fallback for off-grid :math:`\tau`.
A zero prefactor also short-circuits the :math:`K^m_\tau(\xi)` evaluation, which is the common case on the axis.

The :math:`\tau` integral is done with the composite Simpson rule on unit panels :math:`[0,1], [1,2], \ldots, [99,100]`, stopping when a panel contributes less than :math:`10^{-9}` of the running total.
The upper limit of 100 is set by the conical function library.
The conical functions themselves come from the ``Conical`` module of :cite:p:`DUNSTER2017193` and the complex Gamma function from :cite:p:`zhang1996computation`.

``Calculate_F`` and ``Calculate_Sm``, which evaluate the corresponding sum over :math:`S_m` for the potential, are currently unused; the charge density needs only the :math:`H_m` sum.

Changes made to the original program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The program as first written returned NaN at every grid point. The following were corrected.

* **Single precision truncation in the Gamma function.** ``mygamma`` assembled its result with ``CMPLX(GR, GI)`` and no ``kind=``, so the value was truncated to single precision. For :math:`\tau \approx 70` the magnitude :math:`\left| \Gamma\!\left(\frac{1}{2} \pm m + i\tau\right) \right| \sim e^{-\pi\tau/2} \sim 10^{-48}` is perfectly representable in double precision but underflows to zero in single. Numerator and denominator both became zero, so the ratio was :math:`0/0`, and the resulting NaN propagated through the whole integral. This was the source of the NaN.
* **The** :math:`m = 0` **term was skipped.** The summation ran from :math:`m = 1`, although the :math:`\left(2 - \delta_{0m}\right)` factor is only meaningful when :math:`m = 0` is included. For a charge on the axis, where :math:`\xi_c = 1`, azimuthal symmetry means :math:`m = 0` is the *only* term that contributes, so the result was identically zero once the NaN was dealt with.
* **The** :math:`\left(2 - \delta_{0m}\right)` **factor was lost to operator precedence.** It was written as ``(-1.0d0)**m**(2.0d0 - delta_k(0, m))``, and since ``**`` is right-associative in Fortran this parses as :math:`(-1)^{\left(m^{2-\delta_{0m}}\right)}`, which reduces to :math:`(-1)^m` and so halves every :math:`m \geq 1` term. It is also nonconforming for a negative real base with a real exponent. The on-axis charge masked it, since :math:`K^m_\tau(\xi_c) = 0` for :math:`m \geq 1` when :math:`\xi_c = 1`.
* **The integrand was evaluated at** :math:`\tau = 0`. The first Simpson panel starts at :math:`\tau = 0`, but the conical function library requires :math:`\tau > 0`. Since :math:`H_m` carries an explicit :math:`\tau \tanh(\pi\tau)` factor the integrand vanishes there, and the limit is now taken directly instead of calling the library outside its domain.
* **The conical function library has a gap just above** :math:`x = 1`. It returns ``ierr = 1`` for :math:`0 < x - 1 \lesssim 4.6 \times 10^{-4}`, and at :math:`x = 1` exactly it returns zero with ``ierr = 0``, a wrong value with no error flag. This matters because the source charge sits on the axis at :math:`\xi_c = 1` and field points near the axis fall inside the gap. It is covered by ``conicp_series``, described below.
* **A race in the OpenMP grid loop.** ``ierr`` in ``Calculate_Hm`` and ``Calculate_Sm`` was declared with an initializer, which in Fortran implies ``SAVE``, so the variable was shared between threads; the ``ierr == 2`` path could then stop the run spuriously.
* **Smaller fixes.** The electron-position printout read the uninitialized ``phi_p`` instead of ``phi_c``; the interval-convergence test is now guarded against a zero running total; a failed file open stops the program instead of writing to an undefined unit; and the makefile linked ``bessel.o``, ``gammaCHI.o`` and ``air.o`` without listing them as prerequisites, so a partial build could link stale objects.

Conical functions near :math:`x = 1`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For :math:`x - 1 < 10^{-3}` the conical functions are evaluated by ``conicp_series`` rather than by the library, using the Gauss series

.. math::
    P^\mu_{-\frac{1}{2}+i\tau}(x) = \frac{C_\mu}{\Gamma(1+\mu)}
    \left(\frac{x-1}{x+1}\right)^{\frac{\mu}{2}}
    {}_2F_1\!\left(\tfrac{1}{2} - i\tau,\; \tfrac{1}{2} + i\tau;\; 1 + \mu;\; \tfrac{1-x}{2}\right)\, ,
    \qquad
    C_\mu = \prod_{k=0}^{\mu-1}\left[ \left(k + \tfrac{1}{2}\right)^2 + \tau^2 \right]\, .

Both Pochhammer products are real, so this needs no complex arithmetic and contains nothing that can underflow the way the complex Gamma ratio does.
With :math:`|1-x|/2 < 5\times10^{-4}` it converges in at most about 13 terms, and it reproduces the exact limit at :math:`x = 1` on its own: :math:`P^0 = 1` and :math:`P^\mu = 0` for :math:`\mu \geq 1`.

The threshold :math:`x - 1 = 10^{-3}` sits inside the range where both the series and the library are reliable.
Checked against a quadruple precision evaluation of the same series over the range in which it is used (:math:`10^{-8} \leq x - 1 \leq 10^{-3}`, :math:`0 \leq \mu \leq 30`, :math:`\tau \leq 100`), the series is accurate to :math:`3 \times 10^{-15}`, against about :math:`10^{-10}` for the library.
The series is only used near :math:`x = 1` because it loses accuracy to cancellation once :math:`\tau^2 |1-x|` becomes large.

Truncation of the :math:`m` summation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Reference :cite:p:`peridier:4888` limits the sum to :math:`m \leq 30` and truncates it where the terms start to alternate in sign and grow exponentially.
This is applied in ``Truncate_Sum``, which stops the sum when the magnitude of the terms rises on several consecutive steps --- fewer if the sign alternates as well --- and then keeps the terms up to the *smallest* one, where the truncation error of an asymptotic series is smallest.
A sum that simply converges instead exits on a relative tolerance against the largest term.
The earlier version compared successive differences, :math:`|T_{m-1} - T_m| > 10\,|T_{m-2} - T_{m-1}|`, which never looked at the signs, could be tripped by a single noisy term, and compared two vanishingly small differences once the sum had converged.

.. index:: Verlet, Beeman, code, Fowler-Nordheim, Metropolis-Hastings, image charge, conical functions, prolate spheroidal