//
// Created by Jake on 14/02/2021.
//

#ifndef NBODYPROBLEM_QUADTREE_H
#define NBODYPROBLEM_QUADTREE_H

#include <vector>
#include "vector2.h"


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

enum quadrant {
    TL,
    TR,
    BL,
    BR
};

class QuadNode {
public:
    virtual void insert(Particle* p){};
};

class QuadLeaf : QuadNode {
public:
    Particle *p;

    QuadLeaf(){
        p = nullptr;
    }

    void insert(Particle* p){
        this->p = p;
    }
};

class QuadTree : QuadNode {
    Vector2 TL_bound;
    Vector2 BR_bound;
    Vector2 Centroid;

    QuadNode *topLeft;
    QuadNode *topRight;
    QuadNode *botLeft;
    QuadNode *botRight;

public:
    float Mass; // Total mass of children

    QuadTree(Vector2 TL_bound, Vector2 BR_bound){
        this->TL_bound = TL_bound;
        this->BR_bound = BR_bound;
        this->Centroid = (TL_bound + BR_bound) / 2;
        this->Mass = 0;
        this->topLeft = nullptr;
        this->topRight = nullptr;
        this->botLeft = nullptr;
        this->botRight = nullptr;
    }

    QuadTree(std::vector<Particle> particles){
        this->setBounds(particles);
        this->Centroid = (TL_bound + BR_bound) / 2;
        this->Mass = 0;
        this->topLeft = nullptr;
        this->topRight = nullptr;
        this->botLeft = nullptr;
        this->botRight = nullptr;


        for (size_t i = 0; i < particles.size(); i++) {
            this->insert(&particles[i]);
        }
    }

    void setBounds(std::vector<Particle> particles);
    quadrant findQuadrant(Particle* p);
    void insertIntoSubtree(Particle* p);
    void insert(Particle* p);
};

void QuadTree::setBounds(std::vector<Particle> particles){
    this->TL_bound = particles[0].Position;
    this->BR_bound = particles[0].Position;

    for (size_t i = 1; i < particles.size(); i++) {
        if(particles[i].Position.X < TL_bound.X)
            TL_bound.X = particles[i].Position.X;
        if(particles[i].Position.Y > TL_bound.Y)
            TL_bound.Y = particles[i].Position.Y;
        if(particles[i].Position.X > BR_bound.X)
            BR_bound.X = particles[i].Position.X;
        if(particles[i].Position.Y < BR_bound.Y)
            BR_bound.Y = particles[i].Position.Y;
    }
}

quadrant QuadTree::findQuadrant(Particle *p) {
    float p_X = p->Position.X;
    float p_Y = p->Position.Y;

    if(p_X < Centroid.X && p_Y > Centroid.Y)
        return TL;
    else if(p_X > Centroid.X && p_Y > Centroid.Y)
        return TR;
    else if(p_X < Centroid.X && p_Y < Centroid.Y)
        return BL;
    else
        return BR;
}

void QuadTree::insertIntoSubtree(Particle *p) {
    switch (findQuadrant(p)){
        case TL:
            if(topLeft == nullptr)
                topLeft = new QuadNode();
            else
                topLeft = new QuadTree(TL_bound, Centroid);
            topLeft->insert(p);
            break;
        case TR:
            if(topRight == nullptr)
                topRight = new QuadNode();
            else
                topRight = new QuadTree(Vector2(Centroid.X, TL_bound.Y), Vector2(BR_bound.X, Centroid.Y));
            topRight->insert(p);
            break;
        case BL:
            if(botLeft == nullptr)
                botLeft = new QuadNode();
            else
                botLeft = new QuadTree(Vector2(TL_bound.X, Centroid.Y), Vector2(Centroid.X, BR_bound.Y));
            botLeft->insert(p);
            break;
        case BR:
            if(botRight == nullptr)
                botRight = new QuadNode();
            else
                botRight = new QuadTree(Centroid, BR_bound);
            botRight->insert(p);
            break;
    }
}

void QuadTree::insert(Particle *p) {
    if(p == nullptr)
        return;
    else {
        this->Mass += p->Mass;
        insertIntoSubtree(p);
    }
}

#endif //NBODYPROBLEM_QUADTREE_H
