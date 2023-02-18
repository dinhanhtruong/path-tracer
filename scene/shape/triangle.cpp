#include "triangle.h"

#include "util/CS123Common.h"
#include <iostream>
#include <math.h>

using namespace Eigen;

Triangle::Triangle()
{
}

Triangle::Triangle(Vector3f v1, Vector3f v2, Vector3f v3, Vector3f n1, Vector3f n2, Vector3f n3, int index)
    : _v1(v1), _v2(v2), _v3(v3), _n1(n1), _n2(n2), _n3(n3), m_index(index)
{
    _centroid = (_v1 + _v2 + _v3) / 3.f;
    _bbox.setP(_v1);
    _bbox.expandToInclude(_v2);
    _bbox.expandToInclude(_v3);
}

bool Triangle::getIntersection(const Ray &ray, IntersectionInfo *intersection) const
{
    //https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    Vector3f edge1, edge2, h, s, q;
    float a, f, u, v;
    edge1 = _v2 - _v1;
    edge2 = _v3 - _v1;

    h = ray.d.cross(edge2);
    a = edge1.dot(h);

    if(floatEpsEqual(a, 0)) {
        return false;
    }
    f = 1/a;
    s = ray.o - _v1;
    u = f * s.dot(h);
    if(u < 0.f || u > 1.f) {
        return false;
    }
    q = s.cross(edge1);
    v = f * ray.d.dot(q);
    if(v < 0.f || u + v > 1.f) {
        return false;
    }
    float t = f * edge2.dot(q);
    if(t > FLOAT_EPSILON) {
        intersection->t = t;
        intersection->object = this;
        return true;
    } else {
        return false;
    }
}


Vector3f Triangle::getNormal(const IntersectionInfo &I) const
{
    //Calculate Barycentric coordinates to get interpolated normal
    Vector3f p = I.hit;
    return getNormal(p);
}

Vector3f Triangle::getNormal(const Vector3f &p) const{
    Vector3f v0 = _v2 - _v1;
    Vector3f v1 = _v3 - _v1;
    Vector3f v2 = p - _v1;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.f - v - w;

    Vector3f n = v0.cross(v1);
    //If normals weren't loaded from file, calculate them instead (This will be flat shading, not smooth shading)
    Vector3f n1 = floatEpsEqual(_n1.squaredNorm(), 0) ? n : _n1;
    Vector3f n2 = floatEpsEqual(_n2.squaredNorm(), 0) ? n : _n2;
    Vector3f n3 = floatEpsEqual(_n3.squaredNorm(), 0) ? n : _n3;
    Vector3f interpolated_normal = (u * n1 + v * n2 + w * n3).normalized();
    return interpolated_normal;
}

BBox Triangle::getBBox() const
{
    return _bbox;
}

Vector3f Triangle::getCentroid() const
{
    return (_v1 + _v2 + _v3) / 3.f;
}

int Triangle::getIndex() const
{
    return m_index;
}

tinyobj::material_t Triangle::getMaterial() const
{
    return m_material;
}
tinyobj::MaterialType Triangle::getMaterialType() const {
    return m_materialType;
}

void Triangle::setMaterial(const tinyobj::material_t &material)
{
    m_material = material;
    // set material type based on raw material coeffs. Only one of the following should be true
    bool diffuse = !Vector3f(material.diffuse).isZero();
    bool emission = !Vector3f(material.emission).isZero();
    bool specular = !Vector3f(material.specular).isZero();
    bool specularExponent = material.shininess > 0;
    bool refraction = material.ior != 1; // ior = 1 <=> no refraction

    // light source: non-zero emissive
    // DIELECTRIC_REFRACTION: illum = 5
    // glossy specular: illum = 2 AND nonzero specular AND nonzero diffuse
    // ideal diffuse: nonzero diffuse, zero for others
    // ideal specular: nonzero specular exp (Ns) and specular color, zero others


    if (emission) {
        m_materialType = tinyobj::MaterialType::LIGHT_SOURCE;
    } else if (diffuse && !specular) {
        m_materialType = tinyobj::MaterialType::IDEAL_DIFFUSE;
    } else if (material.illum == 2 && specular && diffuse) {
        m_materialType = tinyobj::MaterialType::GLOSSY_SPECULAR;
    } else if (material.illum == 7 && refraction){
        m_materialType = tinyobj::MaterialType::DIELECTRIC_REFRACTION; // refraction + fresnel reflection
    } else if (specular && specularExponent) {
        m_materialType = tinyobj::MaterialType::IDEAL_SPECULAR;
    } else {
        std::cout << "NO ASSIGNED MATERIAL!" << std::endl;
    }

}

/**
 * @brief Triangle::emittedRadiance
 * @param rayOut outward ray from the current emissive source/light
 * @return
 */
Eigen::Vector3f Triangle::emittedRadiance(const Eigen::Vector3f &rayOut) const {
    // temp: constant L_e based on material
    return Eigen::Vector3f(m_material.emission);
}


/**
 * @brief Triangle::brdf
 * @param w_i incoming light direction (must point outward from surface point)
 * @param w_o outgoing direction
 * @param normal at the surface point
 * @return Vec3f of BRDF values
 */
Eigen::Vector3f Triangle::brdf(const Eigen::Vector3f &w_i, const Eigen::Vector3f &w_o, const Eigen::Vector3f &normal) const {
    if (normal.dot(w_i) < 0 && m_materialType != tinyobj::MaterialType::DIELECTRIC_REFRACTION) {
        // edge case: w_i does not point away from the surface (e.g. from importance sampling on glossy BRDFs)
        return Vector3f::Zero();
    }

    switch(m_materialType) {
        case tinyobj::MaterialType::IDEAL_DIFFUSE:
            return Eigen::Vector3f(m_material.diffuse)/M_PI;
    case tinyobj::MaterialType::GLOSSY_SPECULAR:{
        float n = m_material.shininess; // specular exponent
        Eigen::Vector3f normalizingConst = Eigen::Vector3f(m_material.specular) * (n+2)/(2*M_PI);
        return normalizingConst * pow(reflectAboutNormal(normal, w_i).dot(w_o), n);
    }
    case tinyobj::MaterialType::IDEAL_SPECULAR:{
        // hack: set f_r = 1 so that the same radiance formula can be used for mirrors (estimate dirac delta = 1)
        return Eigen::Vector3f::Ones();
    }
    case tinyobj::MaterialType::DIELECTRIC_REFRACTION:{
        return Eigen::Vector3f::Ones();
    }
    default:
        return Vector3f::Zero();
    }
}

/**
 * @brief getTriangleArea computes the area of the given triangle
 * @param tri
 */
float Triangle::getTriangleArea() {
    // triangle area = half area of parallelogram spanned by any two edges of the triangle
    Eigen::Vector3f AB = _v2 - _v1;
    Eigen::Vector3f AC = _v3 - _v1;
    return AB.cross(AC).norm() / 2.f;
}

/**
 * @brief sampleTrianglePoint uniformly samples a point on the given triangle
 * @param tri
 * @return sampled triangle point
 */
Eigen::Vector3f Triangle::sampleTrianglePoint() {
    // using equation (1) from section 4.2 of https://www.cs.princeton.edu/~funk/tog02.pdf
    // idea: take linear combination of vertex positions weighted by functions of two random numbers r1 and r2
    float sqrt_r1 = std::sqrt(generateRandFloat01());
    float r2 = generateRandFloat01();
    return (1-sqrt_r1)*_v1 + (sqrt_r1*(1-r2))*_v2 + (r2*sqrt_r1)*_v3;
}
