
.. _program_listing_file_gravity.h:

Program Listing for File gravity.h
==================================

|exhale_lsh| :ref:`Return to documentation for file <file_gravity.h>` (``gravity.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #ifndef GRAVITY_H
   #define GRAVITY_H
   
   #include "grid.h"
   #include "tree.h"
   
   // Function prototypes
   void calculate_potential(Grid *grid, TreeNode *tree, double opening_angle, int start_index, int end_index);
   
   #endif // GRAVITY_H
