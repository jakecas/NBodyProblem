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
struct Particle
{
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
void ComputeForces(std::vector<Particle> &p_bodies, float p_gravitationalTerm, float p_deltaT)
{
	Vector2 direction,
		force, acceleration;

	float distance;

	for (size_t j = 0; j < p_bodies.size(); ++j)
	{
		Particle &p1 = p_bodies[j];
	
		force = 0.f, acceleration = 0.f;
	
		for (size_t k = 0; k < p_bodies.size(); ++k)
		{
			if (k == j) continue;
		
			Particle &p2 = p_bodies[k];
			
			// Compute direction vector
			direction = p2.Position - p1.Position;
			
			// Limit distance term to avoid singularities
			distance = std::max<float>( 0.5f * (p2.Mass + p1.Mass), fabs(direction.Length()) );
			
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
void MoveBodies(std::vector<Particle> &p_bodies, float p_deltaT)
{
	for (size_t j = 0; j < p_bodies.size(); ++j)
	{
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

int main(int argc, char **argv)
{
    char file[64]; memset(file, '\0', 64);
    bool output = true;
    int particleCount = 1024;
    float gTerm = 1.f;
    int maxIteration = 1000;
	float deltaT = 0.005f;

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
			
	for (int iteration = 0; iteration < maxIteration; ++iteration){
		ComputeForces(bodies, gTerm, deltaT);
		MoveBodies(bodies, deltaT);

		if(output){
		    fileOutput.str(std::string());
		    fileOutput << "out/nbody_" << iteration << ".txt";
		    PersistPositions(fileOutput.str(), bodies);
		}
	}

	return 0;
}