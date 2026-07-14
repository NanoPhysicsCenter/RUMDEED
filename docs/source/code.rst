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
The target rate :math:`a^*` is set per module: 25% for the planar field emission chains, 50% for the planar thermal-field chains, and 35% (the optimum for a two-dimensional random walk) for the tip, cylindrical tip and torus chains.

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

.. index:: Verlet, Beeman, code, Fowler-Nordheim, Metropolis-Hastings