//----------------------------------------------------------------------------------------------
//	Filename:	nbody.cpp
//	Author:		Keith Bugeja
//----------------------------------------------------------------------------------------------
//  CPS3236 assignment for academic year 2017/2018:
//	Sample naive [O(n^2)] implementation for the n-Body problem.
//----------------------------------------------------------------------------------------------

#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#include "sys/time.h"

#include "vector2.h"
#include "cmd.h"

/*
 * Constant definitions for field dimensions, and particle masses
 */
const int fieldWidth = 1000;
const int fieldHalfWidth = fieldWidth >> 1;
const int fieldHeight = 1000;
const int fieldHalfHeight = fieldHeight >> 1;

const float minBodyMass = 2.5f;
const float maxBodyMassVariance = 5.f;

/*
 * Particle structure
 */
struct Particle {
    Vector2 Position;
    Vector2 Velocity;
    float	Mass;

    Particle(void)
            : Position( ((float)rand()) / RAND_MAX * fieldWidth - fieldHalfWidth,
                        ((float)rand()) / RAND_MAX * fieldHeight - fieldHalfHeight)
            , Velocity( 0.f, 0.f )
            , Mass ( ((float)rand()) / RAND_MAX * maxBodyMassVariance + minBodyMass )
    { }

    Particle(float x, float y, float m)
            : Position(x, y)
            , Velocity(0.f, 0.f)
            , Mass(m)
    {}
};

/*
 * Compute forces of particles exerted on one another
 */
void ComputeForces(std::vector<Particle> &p_bodies, float p_gravitationalTerm, float p_deltaT, size_t start, size_t end) {
    Vector2 direction,
            force, acceleration;

    float distance;

    #pragma omp for schedule(static)
    for (size_t j = start; j < end; ++j) {
        Particle &p1 = p_bodies[j];

        force = 0.f, acceleration = 0.f;
        for (size_t k = 0; k < p_bodies.size(); ++k) {
            if (k == j) continue;

            Particle &p2 = p_bodies[k];

            // Compute direction vector
            direction = p2.Position - p1.Position;

            // Limit distance term to avoid singularities
            distance = std::max<float>(0.5f * (p2.Mass + p1.Mass), fabs(direction.Length()));

            // Accumulate force
            force += direction / (distance * distance * distance) * p2.Mass;
        }

        // Compute acceleration for body
        acceleration = force * p_gravitationalTerm;

        // Integrate velocity (m/s)
        p1.Velocity += acceleration * p_deltaT;
    }
}

/*
 * Update particle positions
 */
void MoveBodies(std::vector<Particle> &p_bodies, float p_deltaT, size_t start, size_t end) {
    #pragma omp for schedule(static)
    for (size_t j = start; j < end; ++j){
        p_bodies[j].Position += p_bodies[j].Velocity * p_deltaT;
    }
}

/*
 * Commit particle masses and positions to file in CSV format
 */
void PersistPositions(const std::string &p_strFilename, std::vector<Particle> &p_bodies)
{
    std::cout << "Writing to file: " << p_strFilename << std::endl;

    std::ofstream output(p_strFilename.c_str());

    if (output.is_open()) {
        for (int j = 0; j < p_bodies.size(); j++) {
            output << 	p_bodies[j].Mass << ", " <<
                   p_bodies[j].Position.Element[0] << ", " <<
                   p_bodies[j].Position.Element[1] << std::endl;
        }

        output.close();
    }
    else
        std::cerr << "Unable to persist data to file:" << p_strFilename << std::endl;

}

int main(int argc, char **argv) {
    // Start timing
    struct timeval start, end;
    gettimeofday(&start, NULL);

    int rc = MPI_Init(&argc, &argv);

    if(rc != MPI_SUCCESS)
        std::cerr << "Unable to start MPI environment. RC: " << rc << std::endl;

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
        std::cout << "Running on " << size << " nodes." << std::endl;

    // Creating Particle MPI_Type
    const int elements = 3;
    int          elt_lengths[elements] = {2, 2, 1};
    MPI_Datatype types[elements] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Datatype mpi_particle_type;
    MPI_Aint     offsets[elements];

    offsets[0] = offsetof(Particle, Position);
    offsets[1] = offsetof(Particle, Velocity);
    offsets[2] = offsetof(Particle, Mass);

    MPI_Type_create_struct(elements, elt_lengths, offsets, types, &mpi_particle_type);
    MPI_Type_commit(&mpi_particle_type);


    // Initializing default values
    char file[64]; memset(file, '\0', 64);
    bool output = true;
    int particleCount = 1024;
    int maxIteration = 1000;
    float gTerm = 1.f;
    float deltaT = 0.005f;

    // Read command line arguments
    for (int i = 1; i < argc; i++){
        switch(getArgSwitch(argv[i])){
            case INPUT_FILE:
                strncpy(file, argv[++i], 64);
                break;
            case ENABLE_OUTPUT:
                if(strncmp(argv[++i], "false", 5) == 0)
                    output = false;
                break;
            case NUM_OF_PARTICLES:
                particleCount = atoi(argv[++i]);
                break;
            case G_CONST:
                gTerm = atof(argv[++i]);
                break;
            case ITERATIONS:
                maxIteration = atoi(argv[++i]);
                break;
            case TIME_STEP:
                deltaT = atof(argv[++i]);
                break;
            default:
                printf("Invalid arg found");
                break;
        }
    }

    std::vector<Particle> bodies;

    // Load or generate particles
    if(rank == 0) {
        if (strlen(file) == 0) {
            for (int bodyIndex = 0; bodyIndex < particleCount; ++bodyIndex)
                bodies.push_back(Particle());
        } else {
            std::string line;
            std::ifstream input(file);

            if (!input.is_open())
                throw std::runtime_error("Input file could not be read.");

            while (std::getline(input, line)) {
                std::stringstream ss_line(line);
                std::string m, x, y;
                std::getline(ss_line, m, ',');
                std::getline(ss_line, x, ',');
                std::getline(ss_line, y, ',');
                bodies.push_back(Particle(stof(x), stof(y), stof(m)));
            }
            particleCount = bodies.size();
        }
    }

    if(rank == 0)
        std::cout << "Broadcasting particle count." << std::endl;

    MPI_Bcast(&particleCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate enough space on the other nodes
    if(rank != 0){
        bodies.resize(particleCount);
    }

    // Calculate the number of particles each node will be given
    // (if not divisible, the last node will have some extra or some less)
    int chunksize = particleCount / size;

    int counts[size];
    for(int i = 0; i < size - 1; i++){
        counts[i] = chunksize;
    }
    counts[size-1] = particleCount - ((size - 1) * chunksize);
    int displs[size];
    displs[0] = 0;
    for(int i = 1; i < size; i++){
        displs[i] = displs[i-1] + counts[i-1];
    }


    std::stringstream fileOutput;

    if(rank == 0)
        std::cout << "Broadcasting particle list to start calculations." << std::endl;

    // Send the starting data to everyone
    MPI_Bcast(&bodies[0], particleCount, mpi_particle_type, 0, MPI_COMM_WORLD);
    // Do NBody calculations and output files if flag is set
    for (int iteration = 0; iteration < maxIteration; ++iteration){
        #pragma omp parallel
        {
            #pragma omp master
            if (iteration == 0)
                printf("Using %d threads for node %d\n", omp_get_num_threads(), rank);
            ComputeForces(bodies, gTerm, deltaT, displs[rank], displs[rank] + counts[rank]);
            #pragma omp barrier
            MoveBodies(bodies, deltaT, displs[rank], displs[rank] + counts[rank]);
        }
        // Gather the chunks together at every node
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &bodies[0], counts, displs, mpi_particle_type, MPI_COMM_WORLD);

        if(rank == 0) {

            if (output) {
                fileOutput.str(std::string());
                fileOutput << "out/nbody_" << iteration << ".txt";
                PersistPositions(fileOutput.str(), bodies);
            }
        }
    }


    MPI_Finalize();

    if(rank == 0) {
        // Stop timing
        gettimeofday(&end, NULL);

        long seconds = end.tv_sec - start.tv_sec;
        long useconds = end.tv_usec - start.tv_usec;
        long mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
        std::cout << "Performed computation for " << file << " in: " << mtime << " ms" << std::endl;
    }

    return 0;
}