.. -*- coding: utf-8; mode: rst -*-

.. _anisotropic_solid:

======================================================
Anisotropic Solid Formulation in SEM3D
======================================================

Mandel Notation, Constitutive Law, and Potential Energy

.. contents:: Contents
   :local:
   :depth: 2


Introduction
============

This document describes the implementation of the anisotropic elastic
constitutive law in :program:`SEM3D`, as found in
:file:`SEM3D/SRC/Solid/calcul_forces_solid.inc` (internal force
computation) and :file:`SEM3D/SRC/Solid/dom_solid.F90` (energy output).
The code uses **Mandel (Kelvin) notation** throughout, which differs from
the classical Voigt convention by a :math:`\sqrt{2}` scaling of the shear
components; this choice makes the stiffness tensor a symmetric,
positive-definite matrix under the standard Euclidean inner product.


Notation
========

Index ordering
--------------

Both the strain and stress 6-vectors follow the ordering:

.. math::

   1 \leftrightarrow xx,\quad
   2 \leftrightarrow yy,\quad
   3 \leftrightarrow zz,\quad
   4 \leftrightarrow yz,\quad
   5 \leftrightarrow xz,\quad
   6 \leftrightarrow xy.


Mandel strain and stress vectors
---------------------------------

Given the displacement gradient components
:math:`D_{IJ} = \partial u_I / \partial x_J`, the *Mandel strain vector*
:math:`\tilde{\boldsymbol{\varepsilon}}` is defined as:

.. math::

   \tilde{\boldsymbol{\varepsilon}} =
   \begin{pmatrix}
     \varepsilon_{xx} \\[4pt]
     \varepsilon_{yy} \\[4pt]
     \varepsilon_{zz} \\[4pt]
     \sqrt{2}\,\varepsilon_{yz} \\[4pt]
     \sqrt{2}\,\varepsilon_{xz} \\[4pt]
     \sqrt{2}\,\varepsilon_{xy}
   \end{pmatrix}
   =
   \begin{pmatrix}
     D_{xx} \\[4pt]
     D_{yy} \\[4pt]
     D_{zz} \\[4pt]
     (D_{yz}+D_{zy})/\sqrt{2} \\[4pt]
     (D_{xz}+D_{zx})/\sqrt{2} \\[4pt]
     (D_{xy}+D_{yx})/\sqrt{2}
   \end{pmatrix},

where :math:`\varepsilon_{ij} = \tfrac{1}{2}(\partial u_i/\partial x_j +
\partial u_j / \partial x_i)` are the symmetric strain tensor components.

The *Mandel stress vector* :math:`\tilde{\boldsymbol{\sigma}}` is defined
analogously:

.. math::

   \tilde{\boldsymbol{\sigma}} =
   \begin{pmatrix}
     \sigma_{xx} \\[4pt]
     \sigma_{yy} \\[4pt]
     \sigma_{zz} \\[4pt]
     \sqrt{2}\,\sigma_{yz} \\[4pt]
     \sqrt{2}\,\sigma_{xz} \\[4pt]
     \sqrt{2}\,\sigma_{xy}
   \end{pmatrix}.

**Comparison with Voigt notation.**
The classical Voigt engineering convention uses a factor of 2 (not
:math:`\sqrt{2}`) for the shear strains
(:math:`\gamma_{ij} = 2\varepsilon_{ij}`) and leaves the stresses unscaled.
Mandel notation applies :math:`\sqrt{2}` symmetrically to both stress and
strain, so that the stiffness matrix is fully symmetric and the Euclidean
dot product :math:`\tilde{\boldsymbol{\sigma}} \cdot
\tilde{\boldsymbol{\varepsilon}}` equals twice the elastic energy density
(see :ref:`sec_energy`).


Anisotropic Stiffness Matrix
=============================

.. _sec_mandel_matrix:

Full 6×6 Mandel matrix
-----------------------

The constitutive relation in Mandel notation is:

.. math::

   \tilde{\boldsymbol{\sigma}} = \widetilde{\mathbf{C}}\,\tilde{\boldsymbol{\varepsilon}},

where the :math:`6\times 6` Mandel stiffness matrix
:math:`\widetilde{\mathbf{C}}` is:

.. math::

   \widetilde{\mathbf{C}} =
   \begin{pmatrix}
     c_{xxxx}             & c_{xxyy}             & c_{xxzz}             &
       \sqrt{2}\,c_{xxyz} & \sqrt{2}\,c_{xxxz}   & \sqrt{2}\,c_{xxxy}   \\[4pt]
     c_{yyxx}             & c_{yyyy}             & c_{yyzz}             &
       \sqrt{2}\,c_{yyyz} & \sqrt{2}\,c_{yyxz}   & \sqrt{2}\,c_{yyxy}   \\[4pt]
     c_{zzxx}             & c_{zzyy}             & c_{zzzz}             &
       \sqrt{2}\,c_{zzyz} & \sqrt{2}\,c_{zzxz}   & \sqrt{2}\,c_{zzxy}   \\[4pt]
     \sqrt{2}\,c_{yzxx}   & \sqrt{2}\,c_{yzyy}   & \sqrt{2}\,c_{yzzz}   &
       2\,c_{yzyz}        & 2\,c_{yzxz}          & 2\,c_{yzxy}          \\[4pt]
     \sqrt{2}\,c_{xzxx}   & \sqrt{2}\,c_{xzyy}   & \sqrt{2}\,c_{xzzz}   &
       2\,c_{xzyz}        & 2\,c_{xzxz}          & 2\,c_{xzxy}          \\[4pt]
     \sqrt{2}\,c_{xyxx}   & \sqrt{2}\,c_{xyyy}   & \sqrt{2}\,c_{xyzz}   &
       2\,c_{xyyz}        & 2\,c_{xyxz}          & 2\,c_{xyxy}
   \end{pmatrix}.

This matrix is symmetric
(:math:`\widetilde{C}_{IJ}=\widetilde{C}_{JI}`), positive definite for a
physically admissible elastic material, and identical in structure to the
matrix output by the :program:`homofft` program (see
:file:`homofft/doc/homofft.tex`).


Upper-triangle storage (21-component array)
--------------------------------------------

Because :math:`\widetilde{\mathbf{C}}` is symmetric, only 21 independent
components are stored.  In the code, these are packed into
``dom%Cij_(0:20,...)`` in row-major order of the upper triangle:

.. math::

   \begin{array}{cl}
     \texttt{CC(0)}  & \widetilde{C}_{11} = c_{xxxx}                   \\
     \texttt{CC(1)}  & \widetilde{C}_{12} = c_{xxyy}                   \\
     \texttt{CC(2)}  & \widetilde{C}_{13} = c_{xxzz}                   \\
     \texttt{CC(3)}  & \widetilde{C}_{14} = \sqrt{2}\,c_{xxyz}         \\
     \texttt{CC(4)}  & \widetilde{C}_{15} = \sqrt{2}\,c_{xxxz}         \\
     \texttt{CC(5)}  & \widetilde{C}_{16} = \sqrt{2}\,c_{xxxy}         \\
     \texttt{CC(6)}  & \widetilde{C}_{22} = c_{yyyy}                   \\
     \texttt{CC(7)}  & \widetilde{C}_{23} = c_{yyzz}                   \\
     \texttt{CC(8)}  & \widetilde{C}_{24} = \sqrt{2}\,c_{yyyz}         \\
     \texttt{CC(9)}  & \widetilde{C}_{25} = \sqrt{2}\,c_{yyxz}         \\
     \texttt{CC(10)} & \widetilde{C}_{26} = \sqrt{2}\,c_{yyxy}         \\
     \texttt{CC(11)} & \widetilde{C}_{33} = c_{zzzz}                   \\
     \texttt{CC(12)} & \widetilde{C}_{34} = \sqrt{2}\,c_{zzyz}         \\
     \texttt{CC(13)} & \widetilde{C}_{35} = \sqrt{2}\,c_{zzxz}         \\
     \texttt{CC(14)} & \widetilde{C}_{36} = \sqrt{2}\,c_{zzxy}         \\
     \texttt{CC(15)} & \widetilde{C}_{44} = 2\,c_{yzyz}                \\
     \texttt{CC(16)} & \widetilde{C}_{45} = 2\,c_{yzxz}                \\
     \texttt{CC(17)} & \widetilde{C}_{46} = 2\,c_{yzxy}                \\
     \texttt{CC(18)} & \widetilde{C}_{55} = 2\,c_{xzxz}                \\
     \texttt{CC(19)} & \widetilde{C}_{56} = 2\,c_{xzxy}                \\
     \texttt{CC(20)} & \widetilde{C}_{66} = 2\,c_{xyxy}
   \end{array}


Isotropic limit
---------------

For an isotropic material with Lamé parameters :math:`\lambda` and
:math:`\mu`, all off-diagonal blocks (coupling normal to shear) vanish and
the matrix reduces to:

.. math::

   \widetilde{\mathbf{C}}^{\text{iso}} =
   \begin{pmatrix}
     \lambda+2\mu & \lambda      & \lambda      & 0    & 0    & 0    \\
     \lambda      & \lambda+2\mu & \lambda      & 0    & 0    & 0    \\
     \lambda      & \lambda      & \lambda+2\mu & 0    & 0    & 0    \\
     0            & 0            & 0            & 2\mu & 0    & 0    \\
     0            & 0            & 0            & 0    & 2\mu & 0    \\
     0            & 0            & 0            & 0    & 0    & 2\mu
   \end{pmatrix}.

Note that the shear diagonal entries are :math:`2\mu` (not :math:`\mu`):
this is a consequence of the Mandel :math:`\sqrt{2}` scaling absorbed into
both stress and strain.


Constitutive law in the code
-----------------------------

In :file:`calcul_forces_solid.inc`, the Mandel strain components are first
computed (``M_SQRT1_2`` :math:`= 1/\sqrt{2}`):

.. math::

   \tilde\varepsilon_4 = (D_{yz}+D_{zy}) / \sqrt{2}, \quad
   \tilde\varepsilon_5 = (D_{xz}+D_{zx}) / \sqrt{2}, \quad
   \tilde\varepsilon_6 = (D_{xy}+D_{yx}) / \sqrt{2}.

The Mandel stress vector is then obtained by the matrix-vector product
:math:`\tilde{\boldsymbol{\sigma}} = \widetilde{\mathbf{C}}\,\tilde{\boldsymbol{\varepsilon}}`.
The physical (non-scaled) shear stresses needed for the internal force
assembly are recovered afterwards:

.. math::

   \sigma_{yz} = \tilde\sigma_4 / \sqrt{2}, \quad
   \sigma_{xz} = \tilde\sigma_5 / \sqrt{2}, \quad
   \sigma_{xy} = \tilde\sigma_6 / \sqrt{2}.


.. _sec_energy:

Potential Energy Density
========================

Exact expression
----------------

The elastic potential energy density is:

.. math::

   W = \frac{1}{2}\,\boldsymbol{\sigma} : \boldsymbol{\varepsilon}
     = \frac{1}{2}\sum_{i,j}\sigma_{ij}\varepsilon_{ij}
     = \frac{1}{2}
       \bigl(
         \sigma_{xx}\varepsilon_{xx}
        +\sigma_{yy}\varepsilon_{yy}
        +\sigma_{zz}\varepsilon_{zz}
        +2\sigma_{yz}\varepsilon_{yz}
        +2\sigma_{xz}\varepsilon_{xz}
        +2\sigma_{xy}\varepsilon_{xy}
       \bigr).


Mandel dot-product form
-----------------------

In Mandel notation the factor of 2 on the shear terms is absorbed into the
vector components, so :math:`W` reduces to a simple Euclidean dot product:

.. math::

   \boxed{W = \frac{1}{2}\,\tilde{\boldsymbol{\sigma}} \cdot \tilde{\boldsymbol{\varepsilon}}
            = \frac{1}{2}\sum_{I=1}^{6}\tilde\sigma_I\,\tilde\varepsilon_I.}

This is correct because:

.. math::

   \tilde\sigma_I\,\tilde\varepsilon_I =
   \begin{cases}
     \sigma_{ii}\,\varepsilon_{ii}                                        & I \in \{1,2,3\}, \\
     (\sqrt{2}\,\sigma_{ij})(\sqrt{2}\,\varepsilon_{ij}) = 2\sigma_{ij}\varepsilon_{ij}
                                                                          & I \in \{4,5,6\}.
   \end{cases}


Implementation in the code
---------------------------

The potential energy is computed in :file:`dom_solid.F90` (routines
``get_solid_dom_var`` and ``get_solid_dom_elem_energy``).
The Mandel stress vector is built via the full matrix-vector product and is
kept in scaled form (no division by :math:`\sqrt{2}` before the energy
calculation).  The Mandel strain vector uses the same scaled shear
components ``EYZ``, ``EXZ``, ``EXY``.  The energy is then the half dot
product following the equation above:

.. code-block:: fortran

   ! Mandel stress (shear components NOT divided by sqrt(2))
   sigma(4) = DXX*CC(3) + ... + EYZ*CC(15) + EXZ*CC(16) + EXY*CC(17)
   sigma(5) = DXX*CC(4) + ... + EYZ*CC(16) + EXZ*CC(18) + EXY*CC(19)
   sigma(6) = DXX*CC(5) + ... + EYZ*CC(17) + EXZ*CC(19) + EXY*CC(20)

   ! Mandel strain (shear components = sqrt(2)*eps_ij)
   epsilon(1)=DXX; epsilon(2)=DYY; epsilon(3)=DZZ
   epsilon(4)=EYZ; epsilon(5)=EXZ; epsilon(6)=EXY

   ! Dot product  W = 0.5 * sum_I sigma_tilde(I)*epsilon_tilde(I)
   U = 0
   do I = 1, 6
       U = U + sigma(I)*epsilon(I)
   end do
   P_energy = 0.5d0 * U


Isotropic verification
----------------------

For an isotropic medium under a pure shear deformation
:math:`\varepsilon_{xy} = \varepsilon_0` (all other components zero):

.. math::

   \tilde\varepsilon_6 = \sqrt{2}\,\varepsilon_0, \quad
   \tilde\sigma_6 = \widetilde{C}_{66}\,\tilde\varepsilon_6
                  = 2\mu \cdot \sqrt{2}\,\varepsilon_0,

   W = \tfrac{1}{2}\,\tilde\sigma_6\,\tilde\varepsilon_6
     = \tfrac{1}{2}(2\mu\sqrt{2}\,\varepsilon_0)(\sqrt{2}\,\varepsilon_0)
     = 2\mu\varepsilon_0^2,

which matches the classical result :math:`W = 2\mu\varepsilon_{xy}^2`.


Relation to homofft Output
==========================

The :program:`homofft` homogenization program outputs the effective
stiffness tensor in exactly the same Mandel notation and index ordering
described here (see :file:`homofft/doc/homofft.tex`, section on *User
defined model subroutines*).  Therefore, the ``Cij_`` array in
:program:`SEM3D` can be populated directly from the ``Cstar`` output file
of :program:`homofft` without any rescaling.
