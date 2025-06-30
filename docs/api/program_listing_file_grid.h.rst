
.. _program_listing_file_grid.h:

Program Listing for File grid.h
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_grid.h>` (``grid.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #ifndef GRID_H
   #define GRID_H
   
   #include "vector.h"
   
   // Represents a single cell in the grid
   typedef struct {
       Vector3D pos; // Position of the cell center
       double density; // Density at the cell center
       double potential; // Potential at the cell center
   } Cell;
   
   // Represents the 3D grid
   typedef struct {
       int size; // Number of cells in each dimension (size x size x size)
       double box_size; // Physical size of the simulation box
       Cell *cells; // Array of cells
   } Grid;
   
   // Function prototypes
   Grid* create_grid(int size, double box_size);
   void free_grid(Grid *grid);
   void initialize_density(Grid *grid);
   
   #endif // GRID_H
