//
// Created by Jake on 02/02/2021.
//

#ifndef NBODYPROBLEM_CMD_H
#define NBODYPROBLEM_CMD_H

#include <string.h>

enum cmd_arg_switch {
    INPUT_FILE,
    ENABLE_OUTPUT,
    NUM_OF_PARTICLES,
    G_CONST,
    ITERATIONS,
    TIME_STEP
};

cmd_arg_switch getArgSwitch(const char* arg){
    if(strncmp(arg, "-f", 2) == 0)
        return INPUT_FILE;
    else if(strncmp(arg, "-o", 2) == 0)
        return ENABLE_OUTPUT;
    else if(strncmp(arg, "-b", 2) == 0)
        return NUM_OF_PARTICLES;
    else if(strncmp(arg, "-g", 2) == 0)
        return G_CONST;
    else if(strncmp(arg, "-i", 2) == 0)
        return ITERATIONS;
    else if(strncmp(arg, "-d", 2) == 0)
        return TIME_STEP;
}

#endif //NBODYPROBLEM_CMD_H
