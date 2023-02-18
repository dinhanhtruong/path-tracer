/**
 * @file CS123Common.h
 *
 * Contains data structures and macros commonly used in CS123.
 */
#pragma once
#include "util/tiny_obj_loader.h"
#include <optional>
#ifndef __CS123COMMON_H__
#define __CS123COMMON_H__
#include <random>
#include <math.h>
#include "Eigen/Core"

#include "Eigen/Geometry"

//// glu.h in different location on macs
//#ifdef __APPLE__
//#include <glu.h>
//#else
//#include <GL/glu.h>
//#endif

// from http://en.wikipedia.org/wiki/Assertion_(computing)
#define COMPILE_TIME_ASSERT(pred) switch(0){case 0:case pred:;}

typedef float REAL;

#define IMAGE_WIDTH 512
#define IMAGE_HEIGHT 512
#define NUM_SAMPLES_PER_PIXEL 200
#define NUM_DIRECT_LIGHTING_SAMPLES 10 // the number of shadow rays cast from each non-mirror intersection point
#define PATH_TERMINATION_PROB 0.2 // for Russian roulette random path termination

// EXTRA OPTIONS
#define USE_BRDF_IMPORTANCE_SAMPLING true  // applies to diffuse & phong BRDFs
#define USE_STRATIFIED_SUBPIXEL_SAMPLING true
#define DEPTH_OF_FIELD -1 // set focal distance to positive value to enable DoF. Good value for CornellBox is in [2.8,4.2]
#define DEPTH_OF_FIELD_APERATURE 0.35 // set size of the lens on which rays are scattered for DoF. Higher values = more defocus blur.

#define MIN(a, b) (a) < (b) ? (a) : (b)
#define MAX(a, b) (a) > (b) ? (a) : (b)



// ---------------------
// Common math utilities
// ---------------------

const float FLOAT_EPSILON = 1e-4f;
const double DOUBLE_EPSILON = 1e-8;



inline bool floatEpsEqual(float a, float b) {
    // If the difference between a and b is less than epsilon, they are equal
    return fabs(a - b) < FLOAT_EPSILON;
}

inline bool doubleEpsEqual(double a, double b)
{
    // If the difference between a and b is less than epsilon, they are equal
    return fabs(a - b) < DOUBLE_EPSILON;
}

inline Eigen::Vector4f vec3Tovec4(const Eigen::Vector3f &v, float w) {
    return Eigen::Vector4f(v.x(), v.y(), v.z(), w);
}

/**
 * @brief sphericalToCartesian
 * @param theta angle from the upward axis (+z), in radians. In [0, pi].
 * @param phi angle from the +x axis, in radians. In [0,2*pi]
 * @return unit vector corresponding to (theta, phi) in cartesian coordinates
 */
inline Eigen::Vector3f sphericalToCartesian(float theta, float phi) {
    return Eigen::Vector3f(
                std::sin(theta)*std::cos(phi),
                std::sin(theta)*std::sin(phi),
                std::cos(theta)
            ).normalized();
}


/**
 * @brief reflectAboutNormal
 * @param normal unit normal vec
 * @param w_i incident ray (must point outward from surface point and be in the same space as the normal)
 * @return w_o, the reflection of w_i off the surface
 */
inline Eigen::Vector3f reflectAboutNormal(Eigen::Vector3f normal, Eigen::Vector3f w_i) {
    return 2*(w_i.dot(normal))*normal - w_i;
}

inline float RGBradianceToLuminance(Eigen::Vector3f radiance) {
    return radiance.dot(Eigen::Vector3f(0.2126f, 0.7152f, 0.0722f));
}

inline float schlickApprox(float cosTheta, float iorIn, float iorOut) {
    float R0 = pow((iorIn - iorOut)/(iorIn + iorOut), 2);
    return R0 + (1 - R0)*pow(1 - cosTheta, 5);
}

/**
 * @brief refract
 * @param normal in the same hemisphere as wi.
 * @param wi incident ray direction (pointing outward from the surface point)
 * @param iorIn
 * @param iorOut
 * @return optional refracted ray direction. Will not contain a value if there is total internal reflection
 */
inline std::optional<Eigen::Vector3f> refract(Eigen::Vector3f normal, Eigen::Vector3f wi, float iorIn, float iorOut) {
    // check for total internal reflection
    float iorRatio = (iorIn/iorOut);
    float cosThetaIncident = normal.dot(wi); // wi points away from surface
    float cosSquaredThetaTransmitted = 1 - pow(iorRatio, 2)*(1 - pow(cosThetaIncident, 2));
    if (cosSquaredThetaTransmitted < 0) {
        return std::nullopt;
    }
    // calculate refracted ray angle using Snell's
    float cosThetaTransmitted = sqrt(cosSquaredThetaTransmitted);
    Eigen::Vector3f refractedDir = iorRatio*(-wi) + (iorRatio*cosThetaIncident - cosThetaTransmitted)*normal;
    return std::optional<Eigen::Vector3f>{refractedDir};
}

/**
 * @brief rotationTo
 * @param n NORMALIZED vector to which +z=(0,0,1) is rotated
 * @return
 */
Eigen::Transform<float, 3, Eigen::Affine> rotationTo(Eigen::Vector3f &n);

/**
 * @param normal
 * @param wi incident ray direction at the intersection point in world space. Should point away from the surface.
 * @param tri pointer to the triangle hit at the intersection point
 * @return tuple of (1) a unit direction vec sampled from a hemisphere whose base is perpendicular to the given normal,
 *              and (2) the probability density of the sample
 */
std::tuple<Eigen::Vector3f, float> sampleDirection(Eigen::Vector3f &normal, Eigen::Vector3f &wi, tinyobj::MaterialType materialType, const tinyobj::material_t &mat);

Eigen::Vector2f sampleCirclePoint(float radius);


float generateRandFloat01();




#endif
