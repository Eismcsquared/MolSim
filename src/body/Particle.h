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

public:
  explicit Particle(int type = 0);

  Particle(const Particle &other);

  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg,
      int type = 0);

  virtual ~Particle();

  const std::array<double, 3> &getX() const;

  const std::array<double, 3> &getV() const;

  const std::array<double, 3> &getF() const;

  const std::array<double, 3> &getOldF() const;

  double getM() const;

  int getType() const;


  /////////////////
  void setX(std::array<double, 3> x_arg);
  void setV(std::array<double,3> v_arg);
  void setF(std::array<double,3> f_arg);
  void setOldF(std::array<double,3> oldf_arg);
  void setM(double m_arg);
  void setType(int type);
 ////////////////////////


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
