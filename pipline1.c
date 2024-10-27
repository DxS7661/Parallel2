#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

void parallel_sieve(int n, int rank, int size) {
    int low = (n / size) * rank + 1; // Range start
    int high = (n / size) * (rank + 1); // Range end
    bool *is_prime = malloc((high - low + 1) * sizeof(bool));
    for (int i = 0; i <= high - low; i++) is_prime[i] = true;

    int limit = sqrt(n) + 1;
    for (int prime = 2; prime <= limit; prime++) {
        // Broadcast each prime from root process
        MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int start = max(prime*prime, (low + prime - 1) / prime * prime);
        for (int j = start; j <= high; j += prime) {
            is_prime[j - low] = false;
        }
    }

    // Gather all results at the root
    if (rank == 0) {
        // Print or store results for all primes
    }

    free(is_prime);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size, n = 1000000; // Example n value
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parallel_sieve(n, rank, size);

    MPI_Finalize();
    return 0;
}
