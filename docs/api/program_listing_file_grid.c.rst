
.. _program_listing_file_grid.c:

Program Listing for File grid.c
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_grid.c>` (``grid.c``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #include <stdio.h>
   #include <stdlib.h>
   #include <math.h>
   #include "grid.h"
   
   // Creates and initializes a new grid
   Grid* create_grid(int size, double box_size) {
       Grid *grid = (Grid*)malloc(sizeof(Grid));
       if (!grid) {
           perror("Failed to allocate grid");
           exit(EXIT_FAILURE);
       }
   
       grid->size = size;
       grid->box_size = box_size;
       grid->cells = (Cell*)malloc(size * size * size * sizeof(Cell));
       if (!grid->cells) {
           perror("Failed to allocate cells");
           free(grid);
           exit(EXIT_FAILURE);
       }
   
       double cell_size = box_size / size;
       for (int i = 0; i < size; ++i) {
           for (int j = 0; j < size; ++j) {
               for (int k = 0; k < size; ++k) {
                   int index = i * size * size + j * size + k;
                   grid->cells[index].pos.x = (i + 0.5) * cell_size;
                   grid->cells[index].pos.y = (j + 0.5) * cell_size;
                   grid->cells[index].pos.z = (k + 0.5) * cell_size;
                   grid->cells[index].density = 0.0;
                   grid->cells[index].potential = 0.0;
               }
           }
       }
   
       return grid;
   }
   
   // Frees the memory allocated for the grid
   void free_grid(Grid *grid) {
       if (grid) {
           free(grid->cells);
           free(grid);
       }
   }
   
   // Initializes the density field with a sample distribution
   void initialize_density(Grid *grid) {
       double center_x = grid->box_size / 2.0;
       double center_y = grid->box_size / 2.0;
       double center_z = grid->box_size / 2.0;
   
       // Parameters for two unequal Gaussian peaks
       double peak1_x = center_x - 0.1;
       double peak1_y = center_y - 0.1;
       double peak1_z = center_z;
       double peak1_amplitude = 1.0;
       double peak1_sigma = 0.05;
   
       double peak2_x = center_x + 0.1;
       double peak2_y = center_y + 0.1;
       double peak2_z = center_z;
       double peak2_amplitude = 0.7;
       double peak2_sigma = 0.07;
   
       // Parameters for a background density that decreases towards edges
       double decay_amplitude = 0.1;
       double decay_sigma = grid->box_size / 3.0;
   
       for (int i = 0; i < grid->size; ++i) {
           for (int j = 0; j < grid->size; ++j) {
               for (int k = 0; k < grid->size; ++k) {
                   int index = i * grid->size * grid->size + j * grid->size + k;
                   double x = grid->cells[index].pos.x;
                   double y = grid->cells[index].pos.y;
                   double z = grid->cells[index].pos.z;
   
                   // Calculate squared distances
                   double dist_sq_peak1 = pow(x - peak1_x, 2) + pow(y - peak1_y, 2) + pow(z - peak1_z, 2);
                   double dist_sq_peak2 = pow(x - peak2_x, 2) + pow(y - peak2_y, 2) + pow(z - peak2_z, 2);
                   double dist_sq_center = pow(x - center_x, 2) + pow(y - center_y, 2) + pow(z - center_z, 2);
   
                   // Calculate density contribution from each component
                   double density_peak1 = peak1_amplitude * exp(-dist_sq_peak1 / (2 * peak1_sigma * peak1_sigma));
                   double density_peak2 = peak2_amplitude * exp(-dist_sq_peak2 / (2 * peak2_sigma * peak2_sigma));
                   double density_decay = decay_amplitude * exp(-dist_sq_center / (2 * decay_sigma * decay_sigma));
   
                   grid->cells[index].density = density_peak1 + density_peak2 + density_decay;
   
                   // Ensure density is non-negative
                   if (grid->cells[index].density < 0.0) {
                       grid->cells[index].density = 0.0;
                   }
               }
           }
       }
   }
