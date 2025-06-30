#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> // Include MPI header
#include "grid.h"
#include "tree.h"
#include "gravity.h"

#define GRID_SIZE 64 // Increased resolution
#define BOX_SIZE 2.0 // Bigger box
#define OPENING_ANGLE 0.5

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Grid *grid = NULL;
    if (rank == 0) {
        // Rank 0 creates and initializes the grid
        grid = create_grid(GRID_SIZE, BOX_SIZE);
        initialize_density(grid);
    } else {
        // Other ranks create an empty grid structure to receive data
        grid = create_grid(GRID_SIZE, BOX_SIZE); // Allocate memory for cells
    }

    // Broadcast the density data from rank 0 to all other processes
    MPI_Bcast(grid->cells, GRID_SIZE * GRID_SIZE * GRID_SIZE * sizeof(Cell), MPI_BYTE, 0, MPI_COMM_WORLD);

    // All processes build their own tree (since density data is now consistent)
    TreeNode *tree = build_tree(grid);

    // Determine the range of cells for this process to calculate potential
    int total_cells = GRID_SIZE * GRID_SIZE * GRID_SIZE;
    int cells_per_process = total_cells / size;
    int start_index = rank * cells_per_process;
    int end_index = (rank == size - 1) ? total_cells : (rank + 1) * cells_per_process;

    // Calculate potential for the assigned range of cells
    calculate_potential(grid, tree, OPENING_ANGLE, start_index, end_index);

    // Gather the calculated potentials from all processes to rank 0
    if (rank == 0) {
        // Rank 0 already has its part of the potential calculated
        for (int i = 1; i < size; ++i) {
            int recv_start_index = i * cells_per_process;
            int recv_end_index = (i == size - 1) ? total_cells : (i + 1) * cells_per_process;
            int num_cells_to_recv = recv_end_index - recv_start_index;

            double *temp_potentials = (double*)malloc(num_cells_to_recv * sizeof(double));
            if (!temp_potentials) {
                perror("Failed to allocate temporary potential buffer");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            MPI_Recv(temp_potentials, num_cells_to_recv, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Copy received potentials into the grid
            for (int k = 0; k < num_cells_to_recv; ++k) {
                grid->cells[recv_start_index + k].potential = temp_potentials[k];
            }
            free(temp_potentials);
        }

        // Output the potential results
        FILE *potential_file = fopen("potential.dat", "w");
        if (!potential_file) {
            perror("Failed to open potential output file");
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        for (int i = 0; i < GRID_SIZE; ++i) {
            for (int j = 0; j < GRID_SIZE; ++j) {
                int index = i * GRID_SIZE * GRID_SIZE + j * GRID_SIZE + GRID_SIZE / 2;
                fprintf(potential_file, "%f %f %f\n", grid->cells[index].pos.x, grid->cells[index].pos.y, grid->cells[index].potential);
            }
            fprintf(potential_file, "\n");
        }

        fclose(potential_file);

        // Output the density results
        FILE *density_file = fopen("density.dat", "w");
        if (!density_file) {
            perror("Failed to open density output file");
            MPI_Finalize();
            return EXIT_FAILURE;
        }

        for (int i = 0; i < GRID_SIZE; ++i) {
            for (int j = 0; j < GRID_SIZE; ++j) {
                int index = i * GRID_SIZE * GRID_SIZE + j * GRID_SIZE + GRID_SIZE / 2;
                fprintf(density_file, "%f %f %f\n", grid->cells[index].pos.x, grid->cells[index].pos.y, grid->cells[index].density);
            }
            fprintf(density_file, "\n");
        }

        fclose(density_file);

        printf("Simulation finished. Output written to potential.dat and density.dat\n");
    } else {
        // Other ranks send their calculated potentials to rank 0
        int num_cells_to_send = end_index - start_index;
        double *temp_potentials = (double*)malloc(num_cells_to_send * sizeof(double));
        if (!temp_potentials) {
            perror("Failed to allocate temporary potential buffer");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Copy calculated potentials into the temporary buffer
        for (int k = 0; k < num_cells_to_send; ++k) {
            temp_potentials[k] = grid->cells[start_index + k].potential;
        }

        MPI_Send(temp_potentials, num_cells_to_send, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        free(temp_potentials);
    }

    // Clean up
    free_grid(grid);
    free_tree(tree);

    MPI_Finalize();

    return 0;
}
