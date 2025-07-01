#ifndef GRAVITY_H
#define GRAVITY_H

#include "grid.h"
#include "tree.h"

// Function prototypes
void calculate_potential(Grid *grid, TreeNode *tree, double opening_angle, int start_index, int end_index);

// Helper to get 1D index from 3D local coordinates
static inline int get_local_idx(const Grid *grid, int i, int j, int k) {
    return (i + grid->nghosts) * grid->local_size_y * grid->local_size_z + j * grid->local_size_z + k;
}

#endif // GRAVITY_H
