/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <string>
/**
 * @brief Class that represents a particle. A particle object contains information about its current position, velocity and mass. Moreover, it stores the force acting on it in the current and the previous time step, which is necessary for update via Velocity-Störmer-Verlet.
 */
class Particle {

private:
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

public:
    explicit Particle(int type = 0);

    Particle(const Particle &other);

    Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
        std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
        int type = 0);

    virtual ~Particle();

    /**
     * The getter for the position of a particle.
     * @return the position of the particle.
     */
    const std::array<double, 3> &getX() const;

    /**
     * The getter for the velocity of a particle.
     * @return the velocity of the particle.
     */
    const std::array<double, 3> &getV() const;

    /**
     * The getter for the force acting on a particle.
     * @return the force acting on the particle.
     */
    const std::array<double, 3> &getF() const;

    /**
     * The getter for the force acting on a particle in the previous time step.
     * @return the force acting on the particle in the previous time step.
     */
    const std::array<double, 3> &getOldF() const;

    /**
     * The getter for the mass of a particle.
     * @return the mass of the particle.
     */
    double getM() const;

    /**
     * The getter for the type of a particle.
     * @return the type of the particle.
     */
    int getType() const;

    /**
     * Determine whether the particle is in the domain of simulation.
     * @return True if the particle is in domain.
     */
    bool isInDomain() const;

    /**
     * The setter for the position of a particle.
     * @param x_arg: The new position of the particle.
     */
    void setX(std::array<double, 3> x_arg);
    /**
     * The setter for the velocity of a particle.
     * @param v_arg: The new position of the particle.
     */
    void setV(std::array<double,3> v_arg);
    /**
     * The setter for the force acting on a particle.
     * @param f_arg: The new force acting on the particle.
     */
    void setF(std::array<double,3> f_arg);
    /**
     * The setter for the force acting on a particle in the previous time step.
     * @param oldf_arg: The new position of the particle in the previous time step.
     */
    void setOldF(std::array<double,3> oldf_arg);
    /**
     * The setter for the mass of a particle.
     * @param m_arg: The new mass of the particle.
     */
    void setM(double m_arg);
    /**
     * The setter for the type of a particle.
     * @param type: The new type of the particle.
     */
    void setType(int type);

    /**
     * Remove the particle from the domain.
     */
    void removeFromDomain();

    /**
     * @brief Compare two particles based on their position, velocity, current force and type, under consideration of a
     * relative error up to 1e-5 and an absolute error up to 1e-12.
     * @param other: The Particle to be compared with this.
     * @return
     */
    bool operator==(const Particle &other) const;

    std::string toString() const;

};

std::ostream &operator<<(std::ostream &stream, Particle &p);
