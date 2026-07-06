.. -*- coding: utf-8; mode: rst -*-

.. _anisotropic_fluid:

======================================================
Anisotropic Fluid Formulation in SEM3D
======================================================

Acoustic Wave Propagation with Anisotropy Carried by the DENSITY

.. contents:: Contents
   :local:
   :depth: 2


Introduction
============

This document describes anisotropic-density acoustic materials in
:program:`SEM3D`'s regular fluid domain (:code:`DM_FLUID_CG`,
:file:`SEM3D/SRC/Fluid/`). Anisotropy is a per-model flag (``dom%aniso``) on
that domain, not a separate domain type — the same pattern the solid domain
already uses for isotropic vs. Hooke-anisotropic elastic materials (there is
no dedicated "SolidAniso" domain either).

Unlike a naive "anisotropic bulk modulus" model, this places the
anisotropy in the **density**, not in the bulk modulus. This is the physically
correct effective behaviour obtained when one homogenises a finely layered or
rough acoustic medium: the homogenisation produces an **anisotropic effective
mass matrix** (an anisotropic inverse density), while the effective bulk
modulus stays scalar.

.. note::

   **Reference.** P. Cance and Y. Capdeville, *Validity of the acoustic
   approximation for elastic waves in heterogeneous media*, **Geophysics**
   80(4), T161–T173, 2015 (doi:10.1190/geo2014-0397.1). The paper shows that
   small-scale heterogeneities give rise to a *natural acoustic effective
   anisotropy through an anisotropic effective mass matrix* — i.e. anisotropy
   in the density. This is the SEM3D counterpart of that result and is
   designed to consume the effective model produced by :program:`homofft`
   (acoustic ``compute_effectiveL_r_3d``).

It reuses the regular fluid domain's **velocity-potential** formulation
unchanged, generalising the scalar inverse density
:math:`\mathrm{IDensity}=1/\rho` to a symmetric positive-definite tensor
(``m_IDensTensor``), allocated on that same domain whenever ``dom%aniso`` is
set.


Governing Equation
===================

The unknown is the velocity potential :math:`\varphi`. With scalar bulk
modulus :math:`\kappa` and symmetric inverse-density tensor
:math:`\rho^{-1}_{ij}`:

.. math::

   \frac{1}{\kappa}\,\frac{\partial^2 \varphi}{\partial t^2}
   = \frac{\partial}{\partial x_i}
     \left( \rho^{-1}_{ij}\,\frac{\partial \varphi}{\partial x_j} \right),
   \qquad
   v_i = \rho^{-1}_{ij}\,\frac{\partial \varphi}{\partial x_j},

with Einstein summation (:math:`i,j \in \{x,y,z\}`). The physical particle
velocity :math:`v_i` is the contraction of the inverse-density tensor with the
potential gradient. The pressure is :math:`p = -\dot\varphi`.

In the isotropic limit (:math:`\rho^{-1}_{ij} = \rho^{-1}\delta_{ij}`) this
reduces to the standard acoustic equation of the regular fluid domain,
:math:`\kappa^{-1}\ddot\varphi = \nabla\cdot(\rho^{-1}\nabla\varphi)`, with
wave speed :math:`c = \sqrt{\kappa/\rho}`.

Comparison with the regular fluid
---------------------------------

============================  ==============================  ==============================
Term                          ``dom%aniso = .false.``         ``dom%aniso = .true.``
============================  ==============================  ==============================
Inertial / mass term          :math:`1/\kappa` (scalar)       :math:`1/\kappa` (scalar)
Spatial operator coefficient  :math:`1/\rho` (scalar)         :math:`\rho^{-1}_{ij}` (tensor)
Mass matrix                   :math:`\int J/\kappa`           :math:`\int J/\kappa` (identical)
Particle velocity             :math:`\rho^{-1}\nabla\varphi`  :math:`\rho^{-1}_{ij}\partial_j\varphi`
============================  ==============================  ==============================

The mass matrix is therefore identical to the regular fluid; the only
generalisation is that the scalar :math:`1/\rho` inside the spatial operator
becomes the inverse-density tensor.


Weak Form and Internal Forces
-----------------------------

Multiplying by a test function and integrating by parts over an element gives
the internal-force integrand assembled in
:file:`SEM3D/SRC/Fluid/calcul_forces_fluid_aniso.inc`:

.. math::

   \mathbf{s} = \boldsymbol{\rho}^{-1}\,\nabla\varphi,
   \qquad
   \begin{pmatrix} s_x \\ s_y \\ s_z \end{pmatrix}
   =
   \begin{pmatrix}
     \rho^{-1}_{11} & \rho^{-1}_{12} & \rho^{-1}_{13} \\
     \rho^{-1}_{12} & \rho^{-1}_{22} & \rho^{-1}_{23} \\
     \rho^{-1}_{13} & \rho^{-1}_{23} & \rho^{-1}_{33}
   \end{pmatrix}
   \begin{pmatrix}
     \partial \varphi/\partial x \\
     \partial \varphi/\partial y \\
     \partial \varphi/\partial z
   \end{pmatrix}.

The flux :math:`\mathbf{s}` (which is also the particle velocity) is projected
back to reference coordinates and accumulated into the element force vector via
the standard SEM transpose-differentiation step — exactly as in the regular
fluid kernel, but with the scalar :math:`1/\rho` replaced by the tensor.


The inverse-density tensor
==========================

Storage convention
------------------

The six independent components of :math:`\rho^{-1}_{ij}` are stored in
``m_IDensTensor(0:5,...)``; the scalar bulk modulus :math:`\kappa` is stored in
``m_Lambda``:

.. math::

   \begin{array}{cl}
     \texttt{m\_IDensTensor(0)} & \rho^{-1}_{11} \\
     \texttt{m\_IDensTensor(1)} & \rho^{-1}_{22} \\
     \texttt{m\_IDensTensor(2)} & \rho^{-1}_{33} \\
     \texttt{m\_IDensTensor(3)} & \rho^{-1}_{12} \\
     \texttt{m\_IDensTensor(4)} & \rho^{-1}_{13} \\
     \texttt{m\_IDensTensor(5)} & \rho^{-1}_{23} \\
     \texttt{m\_Lambda}         & \kappa \\
   \end{array}

The tensor must be symmetric positive-definite for the scheme to be stable.

Relation to wave speeds
-----------------------

For an orthorhombic medium (:math:`\rho^{-1}_{12}=\rho^{-1}_{13}=\rho^{-1}_{23}=0`)
the principal-axis wave speeds are:

.. math::

   V_x = \sqrt{\kappa\,\rho^{-1}_{11}}, \quad
   V_y = \sqrt{\kappa\,\rho^{-1}_{22}}, \quad
   V_z = \sqrt{\kappa\,\rho^{-1}_{33}}.

Equivalently, given speeds and a scalar bulk modulus :math:`\kappa`:

.. math::

   \rho^{-1}_{11} = V_x^2/\kappa, \quad
   \rho^{-1}_{22} = V_y^2/\kappa, \quad
   \rho^{-1}_{33} = V_z^2/\kappa.

The CFL phase speed used in :func:`fluid_aniso_Pspeed` is
:math:`\sqrt{\kappa\,\lambda_{\max}(\boldsymbol{\rho}^{-1})}`, bounded above via
the Gershgorin radius of the symmetric tensor.


Material Input Formats
======================

Three input formats are supported, selected via the ``deftype`` keyword in
``material.spec``.

Fluid_Aniso (constant, isotropic-equivalent)
---------------------------------------------

``deftype = Fluid_Aniso; spacedef = constant;`` with ``rho``/``vp``/``vs``.
This is the special case :math:`\rho^{-1}_{ij}=\rho^{-1}\delta_{ij}`,
:math:`\kappa=\rho V_p^2` (mirrors the plain ``Vp/Vs/Rho`` isotropic material
in :file:`define_arrays.F90`) — useful for verifying the tensor kernel reduces
exactly to the scalar one. Example::

  material 0 {
    domain   = fluid;
    deftype  = Fluid_Aniso;
    spacedef = constant;
    rho = 1000.0;
    vp  = 1500.0;
    vs  = 0.0;
  };

See :file:`SEM3D/TESTS/NON-REGR/TEST_0009_cube_fluid_aniso`.

.. warning::

   The values supplied are the **inverse-density tensor** and the **inverse
   bulk modulus**, not a stiffness tensor and a density. For historical I/O
   compatibility the HDF5 group / ``prop_field`` keys retain the legacy names
   ``K11..K23`` and ``Rho``, but their physical content is now
   :math:`\rho^{-1}_{ij}` and :math:`1/\kappa` respectively. The
   interpretation is fixed in
   :func:`init_material_properties_fluid_aniso_from_file`.

Fluid_Aniso (HDF5 file)
-----------------------

``deftype = Fluid_Aniso;`` (constant index 17)

Seven datasets are read from a single HDF5 file:

============  ==============================================  =======================
Group name    Physical content                                Unit
============  ==============================================  =======================
``K11``       :math:`\rho^{-1}_{11}`                          m³/kg
``K22``       :math:`\rho^{-1}_{22}`                          m³/kg
``K33``       :math:`\rho^{-1}_{33}`                          m³/kg
``K12``       :math:`\rho^{-1}_{12}`                          m³/kg
``K13``       :math:`\rho^{-1}_{13}`                          m³/kg
``K23``       :math:`\rho^{-1}_{23}`                          m³/kg
``Rho``       :math:`1/\kappa` (inverse bulk modulus)         Pa⁻¹
============  ==============================================  =======================

Each group exposes the standard SEM material-file attributes ``xMinGlob``,
``xMaxGlob`` (float64[3]) and a ``samples`` dataset of shape ``(nx,ny,nz)``,
:math:`n\geq 2`.

Example ``material.spec``::

  material 0 {
    domain   = fluid;
    deftype  = Fluid_Aniso;
    spacedef = file;
    filename = "mat_fluid_aniso.h5";
  };

Python generator for an orthorhombic medium with :math:`V_x=1500`,
:math:`V_y=1200`, :math:`V_z=900` m/s and bulk modulus
:math:`\kappa=2.25\times10^{9}` Pa:

.. code-block:: python

   import numpy as np, h5py

   KAPPA = 2.25e9                      # scalar bulk modulus [Pa]
   inv_kappa = 1.0/KAPPA               # stored under legacy key "Rho"
   components = {
       "K11": 1500.0**2/KAPPA,         # rho^{-1}_11 = Vx^2/kappa
       "K22": 1200.0**2/KAPPA,         # rho^{-1}_22 = Vy^2/kappa
       "K33":  900.0**2/KAPPA,         # rho^{-1}_33 = Vz^2/kappa
       "K12": 0.0, "K13": 0.0, "K23": 0.0,
       "Rho": inv_kappa,               # 1/kappa  (NOT density)
   }
   with h5py.File("mat_fluid_aniso.h5", "w") as f:
       for name, value in components.items():
           grp = f.create_group(name)
           grp.attrs["xMinGlob"] = [0., 0., 0.]
           grp.attrs["xMaxGlob"] = [500., 500., 500.]
           grp.create_dataset("samples", data=np.full((2,2,2), value))

Cstar_Fluid (binary file from homofft)
---------------------------------------

``deftype = Cstar_Fluid;`` (constant index 18)

The acoustic ``Cstar`` produced by :program:`homofft` (``Nd = 3``) is read
directly. Per :file:`homo/src/cut_cstar3d.f90` (``get_iso_param_acoustic3d``),
homofft writes, per grid point, the effective **inverse-density tensor**
:math:`\rho^{*-1}_{ij}` (the ``LIJ`` matrix) plus the effective **inverse bulk
modulus** :math:`1/\kappa^{*}` (``invkappa``). The on-disk component order is
the upper-triangle row-major of the symmetric tensor followed by the scalar:

.. math::

   \rho^{*-1}_{11},\; \rho^{*-1}_{12},\; \rho^{*-1}_{13},\;
   \rho^{*-1}_{22},\; \rho^{*-1}_{23},\; \rho^{*-1}_{33},\; 1/\kappa^{*}

mapped internally (``build_prop_files.F90``) to the ``prop_field`` order
:math:`(\rho^{-1}_{11},\rho^{-1}_{22},\rho^{-1}_{33},\rho^{-1}_{12},
\rho^{-1}_{13},\rho^{-1}_{23},\,1/\kappa)`. SEM3D stores
``m_IDensTensor = `` :math:`\rho^{*-1}` and
``m_Lambda = `` :math:`\kappa^{*} = 1/(\text{7th value})`.

This is the contract that lets :program:`homofft` and :program:`SEM3D`
"talk": the same ``Cstar`` written by homogenisation is consumed without
reinterpretation. A round-trip check (values written by homofft vs. values
loaded into ``m_IDensTensor`` / ``m_Lambda``) is the recommended verification.

Example ``material.spec``::

  material 0 {
    domain   = fluid;
    deftype  = Cstar_Fluid;
    spacedef = file;
    filename = "mat_cstar_fluid";
  };


Solid–Fluid Coupling
====================

At a solid–fluid interface the usual conditions hold: continuity of normal
displacement and of traction (pressure). In the velocity-potential
formulation the fluid displacement is :math:`u^f_i = \rho^{-1}_{ij}\partial_j\varphi`,
so the normal projection used in the solid→fluid term involves the **normal
contraction of the inverse-density tensor**, :math:`n_i\rho^{-1}_{ij}`, instead
of the scalar :math:`1/\rho` of the isotropic case; the fluid→solid term uses
the pressure :math:`p=-\dot\varphi`. See
:file:`SEM3D/SRC/solid_fluid_coupling.f90`.

.. note::

   Coupling of the anisotropic fluid through a :abbr:`PML (Perfectly Matched
   Layer)` is not covered in this formulation; only the non-PML interface is
   supported.


Implementation Notes
====================

Source files
------------

=================================================  ==============================================
File                                               Purpose
=================================================  ==============================================
:file:`SEM3D/SRC/Fluid/champs_fluid.f90`           Domain type ``domain_fluid``: ``aniso`` flag,
                                                   ``m_IDensTensor`` (:math:`\rho^{-1}`, aniso only),
                                                   ``m_IDensity`` (:math:`1/\rho`, iso), ``m_Lambda``
                                                   (:math:`\kappa`, shared), champs ``Phi,VelPhi,ForcesFl``
:file:`SEM3D/SRC/Fluid/dom_fluid.F90`              Alloc, material init, mass, forces, Newmark,
                                                   velocity (scalar or tensor per ``dom%aniso``),
                                                   Pspeed, energy diagnostics
:file:`SEM3D/SRC/Fluid/calcul_forces_fluid.F90`    Dispatches on ``dom%aniso`` (outside the
                                                   vectorised loop) to the scalar or tensor kernel
:file:`SEM3D/SRC/Fluid/calcul_forces_fluid_aniso.inc`
                                                   Tensor internal-force kernel (vectorised with VCHUNK)
:file:`SEM3D/SRC/build_prop_files.F90`             Property registration / Cstar reading
:file:`SEM3D/SRC/read_input.f90`                   ``read_material_spec``: a ``Fluid_Aniso``/
                                                   ``Cstar_Fluid`` deftype sets ``Tdomain%aniso=.true.``
                                                   (constant or file); no domain override
:file:`SEM3D/SRC/define_arrays.F90`                Material dispatch (constant and file paths);
                                                   ``..._from_file`` sets :math:`\kappa = 1/(1/\kappa)`
:file:`COMMON/constants.F90`                       ``MATDEF_FLUID_ANISO = 17``, ``CSTAR_FLUID = 18``
=================================================  ==============================================

Material-type identifiers
--------------------------

.. code-block:: fortran

   integer, parameter :: MATDEF_FLUID_ANISO = 17  ! HDF5: rho^{-1} groups + 1/kappa
   integer, parameter :: CSTAR_FLUID        = 18  ! homofft binary Cstar (acoustic)

Domain identifier
-----------------

There is no dedicated domain or ``material.input`` character for the
anisotropic-density fluid (unlike, historically, the retired ``'A'`` char).
Materials always declare ``domain = fluid`` (``DM_FLUID_CG``, ``'F'`` in
``material.input``); setting ``deftype`` to one of the values above in
``material.spec`` is what flips ``Tdomain%aniso`` / ``dom%aniso`` and
activates the tensor kernel for that run.


Test Cases
==========

- :file:`SEM3D/TESTS/NON-REGR/TEST_0009_cube_fluid_aniso` — constant material;
  isotropic-equivalent (:math:`\rho^{-1}_{ij}=\rho^{-1}\delta_{ij}`) so the
  result must match the regular fluid domain.
- :file:`SEM3D/TESTS/NON-REGR/TEST_0010_cube_fluid_aniso_h5` — spatially varying
  material from HDF5; truly anisotropic density. The HDF5 file is generated by
  ``gen_mat_h5.py`` in the test directory (regenerate it with the
  inverse-density / ``1/kappa`` convention above).
