//
// Created by Jake on 14/02/2021.
//

#ifndef NBODYPROBLEM_QUADTREE_H
#define NBODYPROBLEM_QUADTREE_H

#include <vector>
#include <string>
#include "vector2.h"

#define THETA 1


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

    void print(std::string tabs){
        std::cout << tabs << "M:"<< Mass << " X: " << Position.X << " Y:" << Position.Y << std::endl;
    }
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
    virtual Particle* getAsParticle(){};
    virtual float getMass(){};
    virtual Vector2 getPosition(){};
    virtual bool isLeaf(){};
    virtual void print(std::string tabs){};
};

class QuadLeaf : public QuadNode {
public:
    Particle *p;

    QuadLeaf(){
        p = nullptr;
    }

    void insert(Particle* p){
        this->p = p;
    }

    Particle* getAsParticle(){
        return p;
    }

    float getMass(){
        return p->Mass;
    }

    Vector2 getPosition(){
        return p->Position;
    }

    bool isLeaf(){
        return true;
    }

    void print(std::string tabs){
        this->p->print(tabs);
    }
};

class QuadTree : public QuadNode {
    Vector2 TL_bound;
    Vector2 BR_bound;
    Vector2 Centroid;

    QuadNode *topLeft;
    QuadNode *topRight;
    QuadNode *botLeft;
    QuadNode *botRight;

    float AggrMass; // Total mass of children
    Vector2 AggrCentroid; // Summation of mass*position of each child node (divide by AggrMass to get weighted centroid)
    Particle *p_repr;

public:
    QuadTree(Vector2 TL_bound, Vector2 BR_bound){
        this->TL_bound = TL_bound;
        this->BR_bound = BR_bound;
        this->Centroid = (TL_bound + BR_bound) / 2;
        this->AggrMass = 0;
        this->AggrCentroid = Vector2(0, 0);
        this->p_repr = nullptr;
        this->topLeft = nullptr;
        this->topRight = nullptr;
        this->botLeft = nullptr;
        this->botRight = nullptr;
    }

    QuadTree(std::vector<Particle> &particles){
        this->setBounds(particles);
        this->Centroid = (TL_bound + BR_bound) / 2;
        this->AggrMass = 0;
        this->p_repr = nullptr;
        this->topLeft = nullptr;
        this->topRight = nullptr;
        this->botLeft = nullptr;
        this->botRight = nullptr;


        for (size_t i = 0; i < particles.size(); i++) {
            if(particles[i].Mass != 0)
                this->insert(&particles[i]);
        }
    }

    float getMass(){
        return AggrMass;
    }

    Vector2 getPosition(){
        return AggrCentroid / AggrMass;
    }

    Particle* getAsParticle(){
        if(p_repr == nullptr) {
            Vector2 pos = getPosition();
            this->p_repr = new Particle(pos.X, pos.Y, getMass());
        }
        return p_repr;
    }

    QuadNode* getTopLeftQuad(){
        return topLeft;
    }

    QuadNode* getTopRightQuad(){
        return topRight;
    }

    QuadNode* getBotLeftQuad(){
        return botLeft;
    }

    QuadNode* getBotRightQuad(){
        return botRight;
    }

    bool isLeaf(){
        return false;
    }

    bool isFarEnough(Particle p);
    void setBounds(std::vector<Particle> particles);
    quadrant findQuadrant(Particle* p);
    void insertIntoSubtree(Particle* p);
    void insert(Particle* p);

    void print(std::string tabs){
        this->getAsParticle()->print(tabs);

        if(topLeft != nullptr)
            this->topLeft->print(tabs+"\t");

        if(topRight != nullptr)
            this->topRight->print(tabs+"\t");

        if(botLeft != nullptr)
            this->botLeft->print(tabs+"\t");

        if(botRight != nullptr)
            this->botRight->print(tabs+"\t");

    }
};

bool QuadTree::isFarEnough(Particle p) {
    // Distance from aggregate centroid to p
    Vector2 direction = this->getPosition() - p.Position;
    float distance = fabs(direction.Length());

    // Radius of this quadrant
    float radius = ((BR_bound.X - TL_bound.X) +(TL_bound.Y - BR_bound.Y))/2.f ;

    // When THETA=1, this is checking that the distance is greater than the radius,
    // i.e. the particle is outside the bounds of this quadrant.
    return (radius / distance) < THETA;
}

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
                topLeft = new QuadLeaf();
            else if(topLeft->isLeaf()){
                Particle *tmp = topLeft->getAsParticle();
                delete topLeft;
                topLeft = new QuadTree(TL_bound, Centroid);
                topLeft->insert(tmp);
            }
            topLeft->insert(p);
            break;
        case TR:
            if(topRight == nullptr)
                topRight = new QuadLeaf();
            else if(topRight->isLeaf()) {
                Particle *tmp = topRight->getAsParticle();
                delete topRight;
                topRight = new QuadTree(Vector2(Centroid.X, TL_bound.Y), Vector2(BR_bound.X, Centroid.Y));
                topRight->insert(tmp);
            }
            topRight->insert(p);
            break;
        case BL:
            if(botLeft == nullptr)
                botLeft = new QuadLeaf();
            else if(botLeft->isLeaf()){
                Particle *tmp = botLeft->getAsParticle();
                delete botLeft;
                botLeft = new QuadTree(Vector2(TL_bound.X, Centroid.Y), Vector2(Centroid.X, BR_bound.Y));
                botLeft->insert(tmp);
            }
            botLeft->insert(p);
            break;
        case BR:
            if(botRight == nullptr)
                botRight = new QuadLeaf();
            else if(botRight->isLeaf()){
                Particle *tmp = botRight->getAsParticle();
                delete botRight;
                botRight = new QuadTree(Centroid, BR_bound);
                botRight->insert(tmp);
            }
            botRight->insert(p);
            break;
    }
}

void QuadTree::insert(Particle *p) {
    if(p == nullptr)
        return;
    else {
        this->AggrMass += p->Mass;
        this->AggrCentroid += p->Mass * p->Position;
        insertIntoSubtree(p);
    }
}

#endif //NBODYPROBLEM_QUADTREE_H
