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
#include <omp.h>

#include "sys/time.h"

#include "vector2.h"
#include "quadtree.h"
#include "cmd.h"


Vector2 ComputeSingleForceOnP(Particle p1, Particle p2){
    // Compute direction vector
    Vector2 direction = p2.Position - p1.Position;

    // Limit distance term to avoid singularities
    float distance = std::max<float>(0.5f * (p2.Mass + p1.Mass), fabs(direction.Length()));

    // Accumulate force
    return direction / (distance * distance * distance) * p2.Mass;
}

Vector2 ComputeAggrForceOnP(Particle p, QuadNode* node) {
    if(node == nullptr)
        return 0.f;
    else if(node->isLeaf()){
        if(p.Position == node->getPosition())
            // node is this particle itself, skip.
            return 0.f;
        else
            return ComputeSingleForceOnP(p, *(node->getAsParticle()));
    }

    QuadTree* tree = dynamic_cast<QuadTree*>(node);
    if(tree->isFarEnough(p))
        return ComputeSingleForceOnP(p, *(tree->getAsParticle()));
    else
        return ComputeAggrForceOnP(p, tree->getTopLeftQuad())
            + ComputeAggrForceOnP(p, tree->getTopRightQuad())
            + ComputeAggrForceOnP(p, tree->getBotLeftQuad())
            + ComputeAggrForceOnP(p, tree->getBotRightQuad());
}

/*
 * Compute forces of particles exerted on one another
 */
void ComputeForces(std::vector<Particle> &p_bodies, QuadTree* tree, float p_gravitationalTerm, float p_deltaT) {
	Vector2 force, acceleration;

    #pragma omp for schedule(static)
    for (size_t j = 0; j < p_bodies.size(); ++j) {
        Particle &p1 = p_bodies[j];

        force = 0.f, acceleration = 0.f;
        force += ComputeAggrForceOnP(p1, tree);

        // Compute acceleration for body
        acceleration = force * p_gravitationalTerm;

        // Integrate velocity (m/s)
        p1.Velocity += acceleration * p_deltaT;
    }
}

/*
 * Update particle positions
 */
void MoveBodies(std::vector<Particle> &p_bodies, float p_deltaT)
{
    #pragma omp for schedule(static)
	for (size_t j = 0; j < p_bodies.size(); ++j){
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
	
	if (output.is_open())
	{	
		for (int j = 0; j < p_bodies.size(); j++)
		{
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

    char file[64]; memset(file, '\0', 64);
    bool output = true;
    int particleCount = 1024;
    float gTerm = 1.f;
    int maxIteration = 1000;
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

	std::stringstream fileOutput;
	std::vector<Particle> bodies;

	// Load or generate particles
	if(strlen(file) == 0){
        for (int bodyIndex = 0; bodyIndex < particleCount; ++bodyIndex)
            bodies.push_back(Particle());
	} else {
        std::string line;
        std::ifstream input(file);

        if(!input.is_open())
            throw std::runtime_error("Input file could not be read.");

	    while(std::getline(input, line)){
	        std::stringstream ss_line(line);
	        std::string m, x, y;
            std::getline(ss_line, m, ',');
	        std::getline(ss_line, x, ',');
	        std::getline(ss_line, y, ',');
	        bodies.push_back(Particle(stof(x), stof(y), stof(m)));
	    }
	}

	std::cout << "Building initial tree..." << std::endl;
	QuadTree *tree = new QuadTree(bodies);
    std::cout << "Tree done. Printing...\n" << std::endl;

	// Do NBody calculations and output files if flag is set
	for (int iteration = 0; iteration < maxIteration; ++iteration){
        #pragma omp parallel
        {
            #pragma omp master
            if(iteration == 0)
                printf("Using %d threads\n", omp_get_num_threads());
            ComputeForces(bodies, tree, gTerm, deltaT);
            #pragma omp barrier
            MoveBodies(bodies, deltaT);
        }
        tree = new QuadTree(bodies);

		if(output){
		    fileOutput.str(std::string());
		    fileOutput << "out/nbody_" << iteration << ".txt";
		    PersistPositions(fileOutput.str(), bodies);
		}
	}

    // Stop timing
    gettimeofday(&end, NULL);

    long seconds = end.tv_sec - start.tv_sec;
    long useconds = end.tv_usec - start.tv_usec;
    long mtime = ((seconds) * 1000 + useconds / 1000.0) + 0.5;
    std::cout << "Performed computation for " << file << " in: " << mtime << " ms" << std::endl;

    delete tree;
	return 0;
}