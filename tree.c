#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tree.h"

// Function to create a new tree node
TreeNode* create_tree_node(Vector3D center, double size) {
    TreeNode *node = (TreeNode*)malloc(sizeof(TreeNode));
    if (!node) {
        perror("Failed to allocate tree node");
        exit(EXIT_FAILURE);
    }

    node->center = center;
    node->size = size;
    node->num_particles = 0;
    node->total_mass = 0.0;
    node->center_of_mass = (Vector3D){0, 0, 0};
    for (int i = 0; i < 8; ++i) {
        node->children[i] = NULL;
    }
    node->particle_indices = NULL;

    return node;
}

// Function to insert a particle into the tree
void insert_particle(TreeNode *node, const Grid *grid, int particle_index) {
    if (node->num_particles == 0) {
        node->particle_indices = (int*)malloc(sizeof(int));
        node->particle_indices[0] = particle_index;
        node->num_particles = 1;
        return;
    }

    if (node->children[0] == NULL) {
        // Create children
        double child_size = node->size / 2.0;
        for (int i = 0; i < 8; ++i) {
            Vector3D child_center = node->center;
            child_center.x += (i & 1) ? child_size / 2.0 : -child_size / 2.0;
            child_center.y += (i & 2) ? child_size / 2.0 : -child_size / 2.0;
            child_center.z += (i & 4) ? child_size / 2.0 : -child_size / 2.0;
            node->children[i] = create_tree_node(child_center, child_size);
        }

        // Move existing particle to a child
        int existing_particle_index = node->particle_indices[0];
        const Cell *existing_cell = &grid->cells[existing_particle_index];
        int existing_child_index = 0;
        if (existing_cell->pos.x > node->center.x) existing_child_index |= 1;
        if (existing_cell->pos.y > node->center.y) existing_child_index |= 2;
        if (existing_cell->pos.z > node->center.z) existing_child_index |= 4;
        insert_particle(node->children[existing_child_index], grid, existing_particle_index);
        free(node->particle_indices);
        node->particle_indices = NULL;
    }

    // Insert new particle into a child
    const Cell *cell = &grid->cells[particle_index];
    int child_index = 0;
    if (cell->pos.x > node->center.x) child_index |= 1;
    if (cell->pos.y > node->center.y) child_index |= 2;
    if (cell->pos.z > node->center.z) child_index |= 4;
    insert_particle(node->children[child_index], grid, particle_index);

    node->num_particles++;
}

// Function to compute mass properties of the tree
void compute_mass_properties(TreeNode *node, const Grid *grid) {
    if (node->num_particles == 1) {
        int particle_index = node->particle_indices[0];
        const Cell *cell = &grid->cells[particle_index];
        node->total_mass = cell->density;
        node->center_of_mass = cell->pos;
        return;
    }

    node->total_mass = 0.0;
    node->center_of_mass = (Vector3D){0, 0, 0};

    for (int i = 0; i < 8; ++i) {
        if (node->children[i]) {
            compute_mass_properties(node->children[i], grid);
            node->total_mass += node->children[i]->total_mass;
            node->center_of_mass.x += node->children[i]->center_of_mass.x * node->children[i]->total_mass;
            node->center_of_mass.y += node->children[i]->center_of_mass.y * node->children[i]->total_mass;
            node->center_of_mass.z += node->children[i]->center_of_mass.z * node->children[i]->total_mass;
        }
    }

    if (node->total_mass > 0) {
        node->center_of_mass.x /= node->total_mass;
        node->center_of_mass.y /= node->total_mass;
        node->center_of_mass.z /= node->total_mass;
    }
}

// Builds the octree from the grid data
TreeNode* build_tree(const Grid *grid) {
    Vector3D root_center = {grid->box_size / 2.0, grid->box_size / 2.0, grid->box_size / 2.0};
    TreeNode *root = create_tree_node(root_center, grid->box_size);

    for (int i = 0; i < grid->size * grid->size * grid->size; ++i) {
        if (grid->cells[i].density > 0) {
            insert_particle(root, grid, i);
        }
    }

    compute_mass_properties(root, grid);

    return root;
}

// Frees the memory allocated for the tree
void free_tree(TreeNode *node) {
    if (!node) return;

    for (int i = 0; i < 8; ++i) {
        free_tree(node->children[i]);
    }

    if (node->particle_indices) {
        free(node->particle_indices);
    }

    free(node);
}
