/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once
#include <vector>
#include <array>
#include <string>
#include "utils/ArrayUtils.h"
/**
 * @brief Class that represents a particle. A particle object contains information about its current position, velocity and mass. Moreover, it stores the force acting on it in the current and the previous time step, which is necessary for update via Velocity-Störmer-Verlet.
 */
class Particle {

private:
    /**
     * The number of particles that have been created.
     */
    inline static int numberParticles = 0;

    /**
     * The identity of the particle, used for identifying neighboring particles in a membrane.
     */
    const int id;

    /**
    * Position of the particle
    */
    std::array<double, 3> x;

    /**
    * Velocity of the particle
    */
    std::array<double, 3> v;

    /**
    * Force effective on this particle
    */
    std::array<double, 3> f;

    /**
    * Force which was effective on this particle
    */
    std::array<double, 3> old_f;

    /**
    * Mass of this particle
    */
    double m;

    /**
    * Type of the particle. Use it for whatever you want (e.g. to separate
    * molecules belonging to different bodies, matters, and so on)
    */
    int type;

    /**
     * Flag that determines whether the particle is in domain of simulation.
     */
    bool inDomain;

    /**
      * Parameter epsilon of the Lennard-Jones potential.
      */
    double epsilon;

    /**
     * Parameter sigma of the Lennard-Jones potential.
     */
    double sigma;

    /**
     * The id of the neighbours of the particle.
     */
    std::vector<int> neighbours;

    /**
     * The id of the diagonal neighbours of the particle.
     */
    std::vector<int> diagonalNeighbours;

    /**
     * Stiffness of the membrane.
     */
    double k;

    /**
     * Average bond length.
     */
    double r0;

    /**
     * Mark a particle as stationary (e.g. for walls)
     */
    const bool stationary;

public:
    explicit Particle(int type = 0);

    Particle(const Particle &other);

    Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
        std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
        int type = 0);

    Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type,
            double epsilon, double sigma, bool stationary = false);

    virtual ~Particle();

    /**
     * The getter for the identity.
     * @return The id of the particle.
     */
    inline const int getId() {return id;}

    /**
     * The getter for the position of a particle.
     * @return the position of the particle.
     */
    inline const std::array<double, 3> &getX() const { return x; };

    /**
     * The getter for the velocity of a particle.
     * @return the velocity of the particle.
     */
    inline const std::array<double, 3> &getV() const { return v; };

    /**
     * The getter for the force acting on a particle.
     * @return the force acting on the particle.
     */
    inline const std::array<double, 3> &getF() const { return f; };

    /**
     * The getter for the force acting on a particle in the previous time step.
     * @return the force acting on the particle in the previous time step.
     */
    inline const std::array<double, 3> &getOldF() const { return old_f; };

    /**
     * The getter for the mass of a particle.
     * @return the mass of the particle.
     */
    inline double getM() const { return m; };

    /**
     * The getter for the type of a particle.
     * @return the type of the particle.
     */
    inline int getType() const { return type; };

    /**
     * Determine whether the particle is in the domain of simulation.
     * @return True if the particle is in domain.
     */
    inline bool isInDomain() const { return inDomain; };

    /**
     * Getter for the parameter epsilon.
     * @return The parameter epsilon.
     */
    inline double getEpsilon() const { return epsilon; };

    /**
     * Getter for the parameter sigma.
     * @return The parameter sigma.
     */
    inline double getSigma() const { return sigma; };

    /**
     * Getter for the stiffness constant.
     * @return The stiffness constant.
     */
    inline double getK() const { return k; };

    /**
     * Getter for the average bond length.
     * @return Te average bond length.
     */
    inline double getR0() const { return r0; };

    /**
     * Determine whether a particle is a direct neighbour of the current particle.
     * @param other The particle that should be tested on neighbourhood.
     * @return True, if the given particle is a direct neighbour.
     */
    inline bool isNeighbour(const Particle &other) {
        return std::find(neighbours.begin(), neighbours.end(), other.id) != neighbours.end();
    }

    /**
     * Determine whether a particle is a diagonal neighbour of the current particle.
     * @param other The particle that should be tested on diagonal neighbourhood.
     * @return True, if the given particle is a diagonal neighbour.
     */
    inline bool isDiagonalNeighbour(const Particle &other) {
        return std::find(diagonalNeighbours.begin(), diagonalNeighbours.end(), other.id) != diagonalNeighbours.end();
    }

    /**
     * The setter for the position of a particle.
     * @param x_arg The new position of the particle.
     */
    inline void setX(std::array<double, 3> x_arg) { this->x = x_arg; };
    /**
     * The setter for the velocity of a particle.
     * @param v_arg The new position of the particle.
     */
    inline void setV(std::array<double,3> v_arg) { this->v = v_arg; };
    /**
     * The setter for the force acting on a particle.
     * @param f_arg The new force acting on the particle.
     */
    inline void setF(std::array<double,3> f_arg) { this->f = f_arg; };
    /**
     * The setter for the force acting on a particle in the previous time step.
     * @param oldf_arg The new position of the particle in the previous time step.
     */
    inline void setOldF(std::array<double,3> oldf_arg) { this->old_f = oldf_arg; };
    /**
     * The setter for the mass of a particle.
     * @param m_arg The new mass of the particle.
     */
    inline void setM(double m_arg) { this->m = m_arg; };
    /**
     * The setter for the type of a particle.
     * @param type The new type of the particle.
     */
    inline void setType(int type_arg) { this->type = type_arg; };

    /**
     * Add a particle to the neighbours of the current particle.
     * @param other The particle to be added.
     */
    inline void addNeighbour(const Particle &other) {
        neighbours.push_back(other.id);
    }

    /**
     * Add a particle to the diagonal neighbours of the current particle.
     * @param other The particle to be added.
     */
    inline void addDiagonalNeighbour(const Particle &other) {
        diagonalNeighbours.push_back(other.id);
    }

    /**
     * Remove the particle from the domain.
     */
    void removeFromDomain();

    /**
     * Compare two particles based on their position, velocity, current force and type, under consideration of a
     * relative error up to 1e-5 and an absolute error up to 1e-12.
     * @param other: The Particle to be compared with this.
     * @return
     */
    bool operator==(const Particle &other) const;

    /**
     * Compare both particles based on their identities.
     * @param other The particle to be compared with.
     * @return True if the both particles have the same identity.
     */
    inline bool is(const Particle &other) { return id == other.id; }

    /**
     * Provide a string representation of the particle including, its position, velocity, force and type.
     * @return A string representation of the particle.
     */
    std::string toString() const;

    /**
     * Update the position of the particle according to Velocity-Störmer-Verlet.
     * @param deltaT The time step of the update.
     */
    inline void updateX(double deltaT) {
        if (inDomain && !stationary) {
            x = x + deltaT * v + (deltaT * deltaT / (2 * m)) * f;
        }
    }

    /**
     * Update the velocity of the particle according to Velocity-Störmer-Verlet.
     * @param deltaT The time step of the update.
     */
    inline void updateV(double deltaT) {
        if (inDomain && !stationary) {
            v = v + (deltaT / (2 * m)) * (f + old_f);
        }
    }

    /**
     * Add force to the particles.
     * @param force The force to be added as an array.
     */
    inline void addForce(std::array<double, 3> force) {
        if (inDomain && !stationary) {
            f = f + force;
        }
    }

};

std::ostream &operator<<(std::ostream &stream, Particle &p);
