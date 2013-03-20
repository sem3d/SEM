
Description des Datasets UNV reconnus




Dataset 2412 ::

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
     94  Thin Shell Linear Quadrilateral
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
     115 Solid Linear Brick
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


Dataset 2411 ::

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


Dataset 780 ::

  Name:   Elements
  Status: Obsolete
  Owner:  Simulation
  Revision Date: 26-SEP-1989
  -----------------------------------------------------------------------
   
  Record 1:        FORMAT(8I10)
                   Field 1       -- element label
                   Field 2       -- fe descriptor id
                   Field 3       -- physical property table bin number
                   Field 4       -- physical property table number
                   Field 5       -- material property table bin number
                   Field 6       -- material property table number
                   Field 7       -- color
                   Field 8       -- number of nodes on element
   
  Record 2:  *** FOR NON-BEAM ELEMENTS ***
                   FORMAT(8I10)
                   Fields 1-n    -- node labels defining element
   
  Record 2:  *** FOR BEAM ELEMENTS ONLY ***
                   FORMAT(3I10)
                   Field 1       -- beam orientation node number
                   Field 2       -- beam fore-end cross section bin number
                   Field 3       -- beam fore-end cross section number
                   Field 4       -- beam  aft-end cross section bin number
                   Field 5       -- beam  aft-end cross section number
   
  Record 3:  *** FOR BEAM ELEMENTS ONLY ***
                   FORMAT(8I10)
                   Fields 1-n    -- node labels defining element
   
  Records 1 and 2 are repeated for each non-beam element in the model.
  Records 1 - 5 are repeated for each beam element in the model.
   
  Example 1:  Solid elements
   
      -1
     780
           1       115         1         1         1         1         8         8
          11        12        13        16        21        20        19        15
           2       113         2         2         1         1         8        16
          31        32        33        34        35        36        37        38
          39        40        41        42        43        44        45        46
           .
           .
           .
         124       115         1         1         1         1         8         8
           9        10        11        15        19        18        17        14
      -1
   
  Example 2:  Beam elements
   
      -1
     780
           1        21         1         1         1         1         7         2
           0         1         1         1         1
           1         2
           2        21         1         1         1         1         7         2
           0         1         1         1         1
           3         4
           3        22         1         3         1         1         7         2
           0         1         1         1         2
           5         6
           .
           .
           .
      -1
   
  ------------------------------------------------------------------------------


Dataset 781 ::

  Name:   Nodes - Double Precision
  Status: Obsolete
  Owner:  Simulation
  Revision Date: 25-MAY-1989
  -----------------------------------------------------------------------
   
  Record 1:        FORMAT(4I10)
                   Field 1       -- node label
                   Field 2       -- definition coordinate system number
                   Field 3       -- displacement coordinate system number
                   Field 4       -- color
  Record 2:        FORMAT(1P3D25.16)
                   Fields 1-3    -- 3-dimensional coordinates of node
                                    in the definition system
   
  Records 1 and 2 are repeated for each node in the model.
   
  Example:
   
      -1
     781
         121         0         0        11
     4.9999998882412910D+00   9.9999997764825821D+00   0.0000000000000000D+00
         122         0         0        11
     5.3124998812563717D+00   9.9999997764825821D+00   0.0000000000000000D+00
         123         0         0        11
     5.6249998742714524D+00   9.9999997764825821D+00   0.0000000000000000D+00
      -1
   
  ----------------------------------------------------------------------------


Dataset 15 ::

  Name:   Nodes
  Status: Obsolete
  Owner:  Simulation
  Revision Date: 30-Aug-1987
  Additional Comments: This dataset is written by I-DEAS Test.
  -----------------------------------------------------------------------
   
               Record 1: FORMAT(4I10,1P3E13.5)
                         Field 1 -    node label
                         Field 2 -    definition coordinate system number
                         Field 3 -    displacement coordinate system number
                         Field 4 -    color
                         Field 5-7 -  3 - Dimensional coordinates of node
                                      in the definition system
   
               NOTE:  Repeat record for each node
   
  ------------------------------------------------------------------------------
