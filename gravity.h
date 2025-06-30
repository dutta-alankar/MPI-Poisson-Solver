#ifndef GRAVITY_H
#define GRAVITY_H

#include "grid.h"
#include "tree.h"

// Function prototypes
void calculate_potential(Grid *grid, TreeNode *tree, double opening_angle, int start_index, int end_index);

#endif // GRAVITY_H
