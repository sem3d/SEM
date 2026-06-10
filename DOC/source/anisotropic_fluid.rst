.. -*- coding: utf-8; mode: rst -*-

.. _anisotropic_fluid:

======================================================
Anisotropic Fluid Formulation in SEM3D
======================================================

Acoustic Wave Equation with Anisotropic Bulk Modulus Tensor

.. contents:: Contents
   :local:
   :depth: 2


Introduction
============

This document describes the implementation of the anisotropic acoustic
domain :code:`DM_FLUID_CG_ANISO` in :program:`SEM3D`, as found in
:file:`SEM3D/SRC/FluidAniso/`.  The domain generalises the standard
acoustic (isotropic fluid) domain by replacing the scalar bulk modulus
:math:`\kappa` with a symmetric positive-definite 3×3 tensor
:math:`\mathbf{K}`.

This formulation is appropriate for media where acoustic wave speed
depends on propagation direction — for example, effective media derived
from heterogeneous microstructures via homogenisation (:program:`homofft`).


Governing Equation
==================

The anisotropic acoustic wave equation for the pressure field :math:`P` is:

.. math::

   \rho \, \frac{\partial^2 P}{\partial t^2}
   = \frac{\partial}{\partial x_i}
     \left( K_{ij} \frac{\partial P}{\partial x_j} \right),

where :math:`\rho` is the density and :math:`\mathbf{K}` is the symmetric
anisotropic bulk modulus tensor.  Einstein summation convention applies
(:math:`i,j \in \{x,y,z\}`).

In the isotropic limit (:math:`K_{ij} = \kappa\,\delta_{ij}`) this reduces
to the standard acoustic wave equation :math:`\rho \ddot P = \kappa\,\nabla^2 P`
with wave speed :math:`c = \sqrt{\kappa/\rho}`.


Weak Form and Internal Forces
------------------------------

The spectral-element discretisation is based on the weak form.  Multiplying
by a test function :math:`\varphi` and integrating by parts over an element
:math:`\Omega_e`:

.. math::

   \int_{\Omega_e}
     \mathbf{K}\,\nabla P \cdot \nabla\varphi \; dV.

The integrand defines the internal force vector.  In the code
(:file:`calcul_forces_fluid_aniso.inc`) this is assembled as:

.. math::

   \mathbf{q} = \mathbf{K}\,\nabla P,
   \qquad
   \text{i.e.}\quad
   \begin{pmatrix} q_x \\ q_y \\ q_z \end{pmatrix}
   =
   \begin{pmatrix}
     K_{11} & K_{12} & K_{13} \\
     K_{12} & K_{22} & K_{23} \\
     K_{13} & K_{23} & K_{33}
   \end{pmatrix}
   \begin{pmatrix}
     \partial P/\partial x \\
     \partial P/\partial y \\
     \partial P/\partial z
   \end{pmatrix}.

The resulting flux :math:`\mathbf{q}` is then projected back to reference
coordinates and accumulated into the element force vector via the standard
SEM transpose-differentiation step.


The Kij Tensor
==============

Storage convention
------------------

The six independent components of :math:`\mathbf{K}` are stored in the
array ``m_Kij(0:5,...)`` following the index ordering:

.. math::

   \begin{array}{cl}
     \texttt{m\_Kij(0)} & K_{11} \\
     \texttt{m\_Kij(1)} & K_{22} \\
     \texttt{m\_Kij(2)} & K_{33} \\
     \texttt{m\_Kij(3)} & K_{12} \\
     \texttt{m\_Kij(4)} & K_{13} \\
     \texttt{m\_Kij(5)} & K_{23} \\
   \end{array}

Relation to wave speeds
-----------------------

For an orthorhombic medium (:math:`K_{12}=K_{13}=K_{23}=0`) the tensor is
diagonal and the wave speeds along the principal axes are:

.. math::

   V_x = \sqrt{\frac{K_{11}}{\rho}}, \quad
   V_y = \sqrt{\frac{K_{22}}{\rho}}, \quad
   V_z = \sqrt{\frac{K_{33}}{\rho}}.

Equivalently, given velocities one sets:

.. math::

   K_{11} = \rho V_x^2, \quad
   K_{22} = \rho V_y^2, \quad
   K_{33} = \rho V_z^2.

Isotropic limit
---------------

When :math:`K_{11}=K_{22}=K_{33}=\kappa` and
:math:`K_{12}=K_{13}=K_{23}=0`, the domain is equivalent to a standard
isotropic fluid with bulk modulus :math:`\kappa`.


Material Input Formats
======================

Two input formats are supported, selected via the ``deftype`` keyword in
``material.spec``.

Fluid_Aniso (HDF5 file)
-----------------------

``deftype = Fluid_Aniso;`` (constant index 17)

All seven properties are read from a single HDF5 file specified by
``filename``.  The file must contain one named group per component:

============  ===================================  =======================
Group name    Description                          Unit
============  ===================================  =======================
``K11``       Diagonal component :math:`K_{11}`   Pa
``K22``       Diagonal component :math:`K_{22}`   Pa
``K33``       Diagonal component :math:`K_{33}`   Pa
``K12``       Off-diagonal :math:`K_{12}`          Pa
``K13``       Off-diagonal :math:`K_{13}`          Pa
``K23``       Off-diagonal :math:`K_{23}`          Pa
``Rho``       Density :math:`\rho`                 kg/m³
============  ===================================  =======================

Each group must expose the standard SEM material-file attributes:

- ``xMinGlob`` — float64 array of size 3, minimum coordinates of the grid
- ``xMaxGlob`` — float64 array of size 3, maximum coordinates of the grid
- ``samples``  — float64 dataset with shape ``(nx, ny, nz)``, :math:`n \geq 2`

Example ``material.spec``::

  material 0 {
    domain   = fluid;
    deftype  = Fluid_Aniso;
    spacedef = file;
    filename = "mat_fluid_aniso.h5";
  };

A minimal Python script to generate the HDF5 file for an orthorhombic
medium with :math:`V_x=1500`, :math:`V_y=1200`, :math:`V_z=900` m/s and
:math:`\rho=1000` kg/m³:

.. code-block:: python

   import numpy as np, h5py

   RHO = 1000.0
   components = {
       "K11": RHO * 1500.0**2,   # 2.25e9 Pa
       "K22": RHO * 1200.0**2,   # 1.44e9 Pa
       "K33": RHO *  900.0**2,   # 8.10e8 Pa
       "K12": 0.0, "K13": 0.0, "K23": 0.0,
       "Rho": RHO,
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

The seven properties are read from a binary Cstar file produced by the
:program:`homofft` homogenisation program, using its acoustic output
(``Nd = 3``).  The binary format is the same as the elastic ``CStar``
format but restricted to the :math:`3\times 3` upper-triangle, giving
7 values per grid point (6 Kij components + density).

The on-disk component order (upper-triangle row-major of the 3×3 tensor)
is:

.. math::

   K_{11},\; K_{12},\; K_{13},\; K_{22},\; K_{23},\; K_{33},\; \rho

which is mapped internally to the ``prop_field`` convention
``(K11, K22, K33, K12, K13, K23, Rho)`` used by the rest of the domain.

Example ``material.spec``::

  material 0 {
    domain   = fluid;
    deftype  = Cstar_Fluid;
    spacedef = file;
    filename = "mat_cstar_fluid.h5";
  };


Implementation Notes
====================

Source files
------------

=================================================  ==============================================
File                                               Purpose
=================================================  ==============================================
:file:`SEM3D/SRC/FluidAniso/dom_fluid_aniso.F90`  Domain type, ``m_Kij``/``m_Rho`` storage,
                                                   ``init_material_properties_fluid_aniso``
:file:`SEM3D/SRC/FluidAniso/calcul_forces_fluid_aniso.inc`
                                                   Internal force kernel (vectorised with VCHUNK)
:file:`SEM3D/SRC/build_prop_files.F90`             Property name registration and file reading
                                                   for ``Fluid_Aniso`` and ``Cstar_Fluid``
:file:`SEM3D/SRC/define_arrays.F90`                Dispatch to ``init_material_properties_fluid_aniso``
                                                   for both ``MATERIAL_CONSTANT`` and
                                                   ``MATERIAL_FILE`` paths
:file:`COMMON/constants.F90`                       Constants ``MATDEF_FLUID_ANISO = 17``,
                                                   ``CSTAR_FLUID = 18``
=================================================  ==============================================

Material-type identifiers
--------------------------

.. code-block:: fortran

   integer, parameter :: MATDEF_FLUID_ANISO = 17  ! HDF5 Kij groups
   integer, parameter :: CSTAR_FLUID        = 18  ! homofft binary Cstar (acoustic)

Domain identifier
-----------------

The character ``'A'`` in ``material.input`` activates the
``DM_FLUID_CG_ANISO`` domain (domain index 7).  When ``material.spec`` is
used instead, set ``domain = fluid`` together with one of the two
``deftype`` values above.


Test Cases
==========

- :file:`SEM3D/TESTS/NON-REGR/TEST_0009_cube_fluid_aniso` — constant
  material defined in ``material.input``; isotropic-equivalent properties.
- :file:`SEM3D/TESTS/NON-REGR/TEST_0010_cube_fluid_aniso_h5` — spatially
  varying material read from an HDF5 file; truly anisotropic with
  :math:`V_x=1500`, :math:`V_y=1200`, :math:`V_z=900` m/s.
  The HDF5 file is generated by ``gen_mat_h5.py`` in the test directory.
