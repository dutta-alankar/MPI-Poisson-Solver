
.. _program_listing_file_tree.h:

Program Listing for File tree.h
===============================

|exhale_lsh| :ref:`Return to documentation for file <file_tree.h>` (``tree.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #ifndef TREE_H
   #define TREE_H
   
   #include "grid.h"
   
   // Represents a node in the octree
   typedef struct TreeNode {
       Vector3D center; // Center of the node
       double size; // Size of the node
       int num_particles;
       double total_mass;
       Vector3D center_of_mass;
       struct TreeNode *children[8]; // Pointers to children nodes
       int *particle_indices; // Indices of particles in the node
   } TreeNode;
   
   // Function prototypes
   TreeNode* build_tree(const Grid *grid);
   void free_tree(TreeNode *node);
   
   #endif // TREE_H
