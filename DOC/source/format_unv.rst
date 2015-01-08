
=============================================
Annexe: Description des Datasets UNV reconnus
=============================================


source: http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse


Chaque bloc d'information (ensemble de données) est délimitée par une chaîne
de séparateur, complètement  sauf la colonne 5 et 6, contenant '' -1 ''.

.. figure:: images/unv_fig001.png
   :scale: 50
   :align: center

Le corps de chaque dataset contient des données qui sont dépendante du
meme dataset. L'enregistrement final de l'ensemble de données contient
une ligne de délimitation contenant ''-1'' dans la colonne 5 et 6.

**Dataset 2411** ::

  Name:   Nodes - Double Precision
  Status: Current
  Owner:  Simulation
  Revision Date: 23-OCT-1992 
  ----------------------------------------------------------------------------
  
  Record 1:        FORMAT(4I10)
                   Field 1       -- node label
                   Field 2       -- export coordinate system number
                   Field 3       -- displacement coordinate system number
                   Field 4       -- color
  Record 2:        FORMAT(1P3D25.16)
                   Fields 1-3    -- node coordinates in the part coordinate
                                    system
   
  Records 1 and 2 are repeated for each node in the model.
   
  Example:
   
      -1
    2411
         121         1         1        11
     5.0000000000000000D+00   1.0000000000000000D+00   0.0000000000000000D+00
         122         1         1        11
     6.0000000000000000D+00   1.0000000000000000D+00   0.0000000000000000D+00
      -1
   
  ----------------------------------------------------------------------------



**Dataset 2412** ::

  Name:   Elements
  Status: Current
  Owner:  Simulation
  Revision Date: 14-AUG-1992
  -----------------------------------------------------------------------
   
  Record 1:        FORMAT(6I10)
                   Field 1       -- element label
                   Field 2       -- fe descriptor id
                   Field 3       -- physical property table number
                   Field 4       -- material property table number
                   Field 5       -- color
                   Field 6       -- number of nodes on element
   
  Record 2:  *** FOR NON-BEAM ELEMENTS ***
                   FORMAT(8I10)
                   Fields 1-n    -- node labels defining element
   
  Record 2:  *** FOR BEAM ELEMENTS ONLY ***
                   FORMAT(3I10)
                   Field 1       -- beam orientation node number
                   Field 2       -- beam fore-end cross section number
                   Field 3       -- beam  aft-end cross section number
   
  Record 3:  *** FOR BEAM ELEMENTS ONLY ***
                   FORMAT(8I10)
                   Fields 1-n    -- node labels defining element
   
  Records 1 and 2 are repeated for each non-beam element in the model.
  Records 1 - 3 are repeated for each beam element in the model.
   
  Example:
   
      -1
    2412
           1        11         1      5380         7         2
           0         1         1
           1         2
           2        21         2      5380         7         2
           0         1         1
           3         4
           3        22         3      5380         7         2
           0         1         2
           5         6
           6        91         6      5380         7         3
          11        18        12
           9        95         6      5380         7         8
          22        25        29        30        31        26        24        23
          14       136         8         0         7         2
          53        54
          36       116        16      5380         7        20
         152       159       168       167       166       158       150       151
         154       170       169       153       157       161       173       172
         171       160       155       156
      -1
  
  FE Descriptor Id definitions
  ____________________________
  
     11  Rod
     21  Linear beam
     22  Tapered beam
     23  Curved beam
     24  Parabolic beam
     31  Straight pipe
     32  Curved pipe
     41  Plane Stress Linear Triangle
     42  Plane Stress Parabolic Triangle
     43  Plane Stress Cubic Triangle
     44  Plane Stress Linear Quadrilateral
     45  Plane Stress Parabolic Quadrilateral
     46  Plane Strain Cubic Quadrilateral
     51  Plane Strain Linear Triangle
     52  Plane Strain Parabolic Triangle
     53  Plane Strain Cubic Triangle
     54  Plane Strain Linear Quadrilateral
     55  Plane Strain Parabolic Quadrilateral
     56  Plane Strain Cubic Quadrilateral
     61  Plate Linear Triangle
     62  Plate Parabolic Triangle
     63  Plate Cubic Triangle
     64  Plate Linear Quadrilateral
     65  Plate Parabolic Quadrilateral
     66  Plate Cubic Quadrilateral
     71  Membrane Linear Quadrilateral
     72  Membrane Parabolic Triangle
     73  Membrane Cubic Triangle
     74  Membrane Linear Triangle
     75  Membrane Parabolic Quadrilateral
     76  Membrane Cubic Quadrilateral
     81  Axisymetric Solid Linear Triangle
     82  Axisymetric Solid Parabolic Triangle
     84  Axisymetric Solid Linear Quadrilateral
     85  Axisymetric Solid Parabolic Quadrilateral
     91  Thin Shell Linear Triangle
     92  Thin Shell Parabolic Triangle
     93  Thin Shell Cubic Triangle
   **94  Thin Shell Linear Quadrilateral (for 2D surface elements)**
     95  Thin Shell Parabolic Quadrilateral
     96  Thin Shell Cubic Quadrilateral
     101 Thick Shell Linear Wedge
     102 Thick Shell Parabolic Wedge
     103 Thick Shell Cubic Wedge
     104 Thick Shell Linear Brick
     105 Thick Shell Parabolic Brick
     106 Thick Shell Cubic Brick
     111 Solid Linear Tetrahedron
     112 Solid Linear Wedge
     113 Solid Parabolic Wedge
     114 Solid Cubic Wedge
   **115 Solid Linear Brick (for 3D solid elements)**
     116 Solid Parabolic Brick
     117 Solid Cubic Brick
     118 Solid Parabolic Tetrahedron
     121 Rigid Bar
     122 Rigid Element
     136 Node To Node Translational Spring
     137 Node To Node Rotational Spring
     138 Node To Ground Translational Spring
     139 Node To Ground Rotational Spring
     141 Node To Node Damper
     142 Node To Gound Damper
     151 Node To Node Gap
     152 Node To Ground Gap
     161 Lumped Mass
     171 Axisymetric Linear Shell
     172 Axisymetric Parabolic Shell
     181 Constraint
     191 Plastic Cold Runner
     192 Plastic Hot Runner
     193 Plastic Water Line
     194 Plastic Fountain
     195 Plastic Baffle
     196 Plastic Rod Heater
     201 Linear node-to-node interface
     202 Linear edge-to-edge interface
     203 Parabolic edge-to-edge interface
     204 Linear face-to-face interface
     208 Parabolic face-to-face interface
     212 Linear axisymmetric interface
     213 Parabolic axisymmetric interface
     221 Linear rigid surface
     222 Parabolic rigin surface
     231 Axisymetric linear rigid surface
     232 Axisymentric parabolic rigid surface
  
  ------------------------------------------------------------------------------


**Dataset 2477** (According to IDEAS docs, dataset 2467 is obsolete and
 is replaced by dataset 2477) ::

  Name:   Permanent Groups
  Status: Current
  Owner:  Simulation
  Revision Date: 24-April-2002
  -----------------------------------------------------------------------
  
  Record 1:        FORMAT(8I10)
                   Field 1       -- group number
                   Field 2       -- active constraint set no. for group
                   Field 3       -- active restraint set no. for group
                   Field 4       -- active load set no. for group
                   Field 5       -- active dof set no. for group
                   Field 6       -- active temperature set no. for group
                   Field 7       -- active contact set no. for group
                 **Field 8       -- number of entities in group**
  
  Record 2:        FORMAT(20A2)
                   Field 1       -- group name **(PhysicalVolume - PhysicalSurface)**
  
  Record 3-N:      FORMAT(8I10)
                   Field 1       -- entity type code
                   Field 2       -- entity tag
                   Field 3       -- entity node leaf id.
                   Field 4       -- entity component/ ham id.
                   Field 5       -- entity type code
                   Field 6       -- entity tag
                   Field 7       -- entity node leaf id.
                   Field 8       -- entity component/ ham id.
  
  Repeat record 3 for all entities as defined by record 1, field 8.
  Records 1 thru n are repeated for each group in the model.
  Entity node leaf id. and the component/ ham id. are zero for all
  entities except "reference point", "reference point series"
  and "coordinate system".
  
            Permanent group entity type codes
  
      Entity Type Code        Entity Description
  
             1                coordinate system
             2                data surface thickness
             3                force on point
             4                force on edge
             5                traction on face
             6                pressure on face
             7                nodes
             8                finite elements
             9                dof sets, dof entities
            10                constraint sets, coupled dofs
            11                constraint sets, mpc equations
            12                restraint sets, nodal displacements
            13                restraint sets, nodal temperatures
            14                load sets, nodal forces
            15                load sets, nodal temperatures
            16                load sets, nodal heat sources/sinks
            17                load sets, face pressures
            18                load sets, edge pressures
            19                load sets, face heat fluxes
            20                load sets, edge heat fluxes
            21                load sets, face heat convections
            22                load sets, edge heat convections
            23                load sets, face heat radiations
            24                load sets, edge heat radiations
            25                load sets, element heat generations
            26                load sets, beam temperatures
            27                trace lines
            28                beam force
            29                beam distributed load
            30                data surface
            31                data curve
            32                displacement on point (restraint)
            33                displacement on edge (restraint)
            34                displacement on surface (restraint)
            35                temperature on point (restraint)
            36                temperature on edge (restraint) 
            37                temperature on face (restraint)
            38                temperature on point (temperature)
            39                temperature on edge (temperature)
            40                temperature on face (temperature)
            41                heat source on point
            42                heat flux on edge
            43                convection on edge
            44                radiation on edge
            45                heat flux on face
            46                convection on face
            47                radiation on face
            48                geometry contact region
            49                fe contact region
            50                contact pair
            51                kinematic dof on point
            52                kinematic dof on edge
            53                kinematic dof on face
            54                element definition
            55                anchor node
            56                edge dependancy mesh definition
            57                fem point connector
            58                fem area connector
            59                vertex
            60                edge
            61                face
            62                region
            63                wireframe connector
            64                wireframe curve
            65                wireframe section
            66                wireframe region
            67                reference point
            68                reference point series
            69                centerpoint
  
  Example:
  
    2477
      -1
      -1
           0         0         0         0         0         0         0         1
  PhysicalSurface0
           8        33         0         0
           1         0         0         0         0         0         0         1
  PhysicalSurface1
           8        38         0         0
           2         0         0         0         0         0         0         1
  PhysicalSurface2
           8        38         0         0
           3         0         0         0         0         0         0         1
  PhysicalSurface3
           8        43         0         0
      -1
           0         0         0         0         0         0         0         1
  PhysicalVolume0
           8        44         0         0
           1         0         0         0         0         0         0         1
  PhysicalVolume1
           8        45         0         0
  -1
  
  -----------------------------------------------------------------------

