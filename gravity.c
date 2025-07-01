#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gravity.h"

#define G 1.0 // Gravitational constant

// Calculates the potential at a single point due to a tree node
double calculate_potential_recursive(Vector3D pos, TreeNode *node, double opening_angle, double box_size) {
    double dx = pos.x - node->center_of_mass.x;
    double dy = pos.y - node->center_of_mass.y;
    double dz = pos.z - node->center_of_mass.z;

    double r = sqrt(dx*dx + dy*dy + dz*dz);

    if (r == 0.0) return 0.0;

    if (node->size / r < opening_angle || node->num_particles == 1) {
        return -G * node->total_mass / r;
    } else {
        double potential = 0.0;
        for (int i = 0; i < 8; ++i) {
            if (node->children[i]) {
                potential += calculate_potential_recursive(pos, node->children[i], opening_angle, box_size);
            }
        }
        return potential;
    }
}

// Calculates the potential for a range of cells in the grid
void calculate_potential(Grid *grid, TreeNode *tree, double opening_angle, int start_index, int end_index) {
    for (int i = 0; i < grid->local_size_x; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int index = get_local_idx(grid, i, j, k);
                grid->cells[index].potential = calculate_potential_recursive(grid->cells[index].pos, tree, opening_angle, grid->box_size);
            }
        }
    }
}
