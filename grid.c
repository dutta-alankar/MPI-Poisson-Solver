#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stddef.h> // For offsetof
#include "grid.h"

// Creates and initializes a new distributed grid
Grid* create_distributed_grid(int global_size, double box_size, int nghosts, MPI_Comm comm) {
    Grid *grid = (Grid*)malloc(sizeof(Grid));
    if (!grid) {
        perror("Failed to allocate grid");
        MPI_Abort(comm, EXIT_FAILURE);
    }

    grid->global_size = global_size;
    grid->box_size = box_size;
    grid->nghosts = nghosts;
    grid->comm = comm;

    MPI_Comm_rank(comm, &grid->rank);
    MPI_Comm_size(comm, &grid->size);

    // Decompose the grid along the x-axis
    int base_size = global_size / grid->size;
    int remainder = global_size % grid->size;

    grid->local_size_x = base_size + (grid->rank < remainder ? 1 : 0);
    grid->offset_x = grid->rank * base_size + (grid->rank < remainder ? grid->rank : remainder);

    grid->local_size_y = global_size;
    grid->local_size_z = global_size;
    grid->offset_y = 0;
    grid->offset_z = 0;

    int total_local_x = grid->local_size_x + 2 * nghosts;
    grid->cells = (Cell*)malloc(total_local_x * grid->local_size_y * grid->local_size_z * sizeof(Cell));
    if (!grid->cells) {
        perror("Failed to allocate cells");
        free(grid);
        MPI_Abort(comm, EXIT_FAILURE);
    }

    double cell_size = box_size / global_size;
    for (int i = 0; i < grid->local_size_x; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int global_i = grid->offset_x + i;
                int local_idx_with_ghosts = get_idx(grid, i + nghosts, j, k);

                grid->cells[local_idx_with_ghosts].pos.x = (global_i + 0.5) * cell_size;
                grid->cells[local_idx_with_ghosts].pos.y = (j + 0.5) * cell_size;
                grid->cells[local_idx_with_ghosts].pos.z = (k + 0.5) * cell_size;
                grid->cells[local_idx_with_ghosts].density = 0.0;
                grid->cells[local_idx_with_ghosts].potential = 0.0;
            }
        }
    }

    return grid;
}

void free_grid(Grid *grid) {
    if (grid) {
        free(grid->cells);
        free(grid);
    }
}

void initialize_density(Grid *grid) {
    double center_x = grid->box_size / 2.0;
    double center_y = grid->box_size / 2.0;
    double center_z = grid->box_size / 2.0;

    double peak1_x = center_x - 0.1, peak1_y = center_y - 0.1, peak1_z = center_z;
    double peak1_amplitude = 1.0, peak1_sigma = 0.05;

    double peak2_x = center_x + 0.1, peak2_y = center_y + 0.1, peak2_z = center_z;
    double peak2_amplitude = 0.7, peak2_sigma = 0.07;

    double decay_amplitude = 0.1, decay_sigma = grid->box_size / 3.0;

    for (int i = 0; i < grid->local_size_x; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int index = get_idx(grid, i + grid->nghosts, j, k);
                double x = grid->cells[index].pos.x;
                double y = grid->cells[index].pos.y;
                double z = grid->cells[index].pos.z;

                double dist_sq_peak1 = pow(x - peak1_x, 2) + pow(y - peak1_y, 2) + pow(z - peak1_z, 2);
                double dist_sq_peak2 = pow(x - peak2_x, 2) + pow(y - peak2_y, 2) + pow(z - peak2_z, 2);
                double dist_sq_center = pow(x - center_x, 2) + pow(y - center_y, 2) + pow(z - center_z, 2);

                double density_peak1 = peak1_amplitude * exp(-dist_sq_peak1 / (2 * peak1_sigma * peak1_sigma));
                double density_peak2 = peak2_amplitude * exp(-dist_sq_peak2 / (2 * peak2_sigma * peak2_sigma));
                double density_decay = decay_amplitude * exp(-dist_sq_center / (2 * decay_sigma * decay_sigma));

                grid->cells[index].density = density_peak1 + density_peak2 + density_decay;
                if (grid->cells[index].density < 0.0) {
                    grid->cells[index].density = 0.0;
                }
            }
        }
    }
}

void exchange_ghost_cells(Grid *grid) {
    if (grid->nghosts == 0) return;

    int rank = grid->rank;
    int size = grid->size;
    MPI_Comm comm = grid->comm;

    int left_neighbor = (rank - 1 + size) % size;
    int right_neighbor = (rank + 1) % size;

    int slab_size = grid->nghosts * grid->local_size_y * grid->local_size_z;
    Cell *send_buffer_left = (Cell*)malloc(slab_size * sizeof(Cell));
    Cell *send_buffer_right = (Cell*)malloc(slab_size * sizeof(Cell));
    Cell *recv_buffer_left = (Cell*)malloc(slab_size * sizeof(Cell));
    Cell *recv_buffer_right = (Cell*)malloc(slab_size * sizeof(Cell));

    // Create a custom MPI datatype for the Cell struct
    MPI_Datatype MPI_Cell;
    int blocklengths[3] = {3, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    displacements[0] = offsetof(Cell, pos);
    displacements[1] = offsetof(Cell, density);
    displacements[2] = offsetof(Cell, potential);
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_Cell);
    MPI_Type_commit(&MPI_Cell);

    // Packing logic can be optimized, but this is clear
    // Pack left-side local cells to send to left neighbor
    for (int i = 0; i < grid->nghosts; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int src_idx = get_idx(grid, i + grid->nghosts, j, k);
                int buf_idx = i * grid->local_size_y * grid->local_size_z + j * grid->local_size_z + k;
                send_buffer_left[buf_idx] = grid->cells[src_idx];
            }
        }
    }

    // Pack right-side local cells to send to right neighbor
    for (int i = 0; i < grid->nghosts; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int src_idx = get_idx(grid, grid->local_size_x + i, j, k);
                int buf_idx = i * grid->local_size_y * grid->local_size_z + j * grid->local_size_z + k;
                send_buffer_right[buf_idx] = grid->cells[src_idx];
            }
        }
    }

    MPI_Status status;
    MPI_Sendrecv(send_buffer_right, slab_size, MPI_Cell, right_neighbor, 0,
                 recv_buffer_left, slab_size, MPI_Cell, left_neighbor, 0,
                 comm, &status);

    MPI_Sendrecv(send_buffer_left, slab_size, MPI_Cell, left_neighbor, 1,
                 recv_buffer_right, slab_size, MPI_Cell, right_neighbor, 1,
                 comm, &status);

    // Unpack received data into ghost cells
    for (int i = 0; i < grid->nghosts; ++i) {
        for (int j = 0; j < grid->local_size_y; ++j) {
            for (int k = 0; k < grid->local_size_z; ++k) {
                int buf_idx = i * grid->local_size_y * grid->local_size_z + j * grid->local_size_z + k;
                // From left neighbor into left ghosts
                grid->cells[get_idx(grid, i, j, k)] = recv_buffer_left[buf_idx];
                // From right neighbor into right ghosts
                grid->cells[get_idx(grid, grid->local_size_x + grid->nghosts + i, j, k)] = recv_buffer_right[buf_idx];
            }
        }
    }

    free(send_buffer_left);
    free(send_buffer_right);
    free(recv_buffer_left);
    free(recv_buffer_right);
    MPI_Type_free(&MPI_Cell);
}
