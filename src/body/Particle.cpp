#include "Particle.h"

#include "utils/ArrayUtils.h"
#include "spdlog/spdlog.h"

Particle::Particle(int type_arg) {
  type = type_arg;
  spdlog::trace("Particle generated!");
  f = {0., 0., 0.};
  old_f = {0., 0., 0.};
  inDomain = true;
}

Particle::Particle(const Particle &other) {
  x = other.x;
  v = other.v;
  f = other.f;
  old_f = other.old_f;
  m = other.m;
  type = other.type;
  inDomain = other.inDomain;
  epsilon = other.epsilon;
  sigma = other.sigma;
  spdlog::trace("Particle generated by copy!");
}

Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
                   double m_arg, int type_arg) : x(x_arg), v(v_arg), f{0., 0., 0.}, old_f{0., 0., 0.}, m(m_arg), type(type_arg), inDomain(
        true), epsilon(5), sigma(1) {
  spdlog::trace("Particle generated!");
}



Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type, double epsilon,
                   double sigma): x(x_arg), v(v_arg), f{0., 0., 0.}, old_f{0., 0., 0.}, m(m_arg), type(type), inDomain(
        true), epsilon(epsilon), sigma(sigma){

}

Particle::~Particle() { spdlog::trace("Particle destructed!"); }

const std::array<double, 3> &Particle::getX() const { return x; }

const std::array<double, 3> &Particle::getV() const { return v; }

const std::array<double, 3> &Particle::getF() const { return f; }

const std::array<double, 3> &Particle::getOldF() const { return old_f; }

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

bool Particle::isInDomain() const {return inDomain;}

double Particle::getEpsilon() const {return epsilon;}

double Particle::getSigma() const {return sigma;}


///////////////////////////////////////////////////
void Particle::setX(std::array<double, 3> x_arg) {x = x_arg; }
void Particle::setV(std::array<double, 3> v_arg) {v = v_arg; }
void Particle::setF(std::array<double, 3> f_arg) {f = f_arg; }
void Particle::setOldF(std::array<double, 3> oldf_arg) {old_f = oldf_arg; }
void Particle::setM(double marg) { m = marg; }
void Particle::setType(int type_arg) { type = type_arg; }
///////////////////////////////////////////////////

void Particle::removeFromDomain() {
    inDomain = false;
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: x:" << x << " v: " << v << " f: " << f
         << " type: " << type;
  return stream.str();
}


bool Particle::operator==(const Particle &other) const {
    double delta_abs = 1e-12;
    double delta_rel = 1e-5;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(x[i] - other.x[i]) > delta_rel * std::max(std::abs(x[i]), std::abs(other.x[i])) + delta_abs) {
            return false;
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (std::abs(v[i] - other.v[i]) > delta_rel * std::max(std::abs(v[i]), std::abs(other.v[i])) + delta_abs) {
            return false;
        }
    }
    for (int i = 0; i < 3; ++i) {
        if (std::abs(f[i] - other.f[i]) > delta_rel * std::max(std::abs(f[i]), std::abs(other.f[i])) + delta_abs) {
            return false;
        }
    }
    if (std::abs(m - other.m) > delta_rel * std::max(std::abs(m), std::abs(other.m)) + delta_abs) {
        return false;
    }

    if (std::abs(epsilon - other.epsilon) > delta_rel * std::max(std::abs(epsilon), std::abs(other.epsilon)) + delta_abs) {
        return false;
    }

    if (std::abs(sigma - other.sigma) > delta_rel * std::max(std::abs(sigma), std::abs(other.sigma)) + delta_abs) {
        return false;
    }
    return type == other.type;
}


std::ostream &operator<<(std::ostream &stream, Particle &p) {
  stream << p.toString();
  return stream;
}
