#include "util/CS123Common.h"
#include <iostream>

std::default_random_engine generator;
std::uniform_real_distribution<float> uniformRealDist(0.0,1.0);


/**
 * @brief rotationTo returns a rotation matrix representing a rotation from oldUp to newUp.
 *          The axis of rotation is perpendicular to the spane planned by oldUp and newUp
 * @param oldUp NORMALIZED vector representing a source reference
 * @param newUp NORMALIZED vector to which oldUp is rotated
 * @return Eigen::Transform encoding a pure rotation from oldUp to newUp
 */
Eigen::Transform<float, 3, Eigen::Affine> rotationTo(Eigen::Vector3f &oldUp, Eigen::Vector3f &newUp) {
    // initialize transformation as identity
    Eigen::Transform<float, 3, Eigen::Affine> t = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    // axis of rotation is normal to both the oldUp and newUp
    Eigen::Vector3f rotAxis = oldUp.cross(newUp).normalized();
    // amount of rotation = angle between +z and n. Always between 0 and pi
    float rotAngle = acos(newUp.dot(oldUp)); // n dot up = cos(theta)
    t.rotate(Eigen::AngleAxisf(rotAngle, rotAxis));

    return t;
}

/**
 * @brief uniformSampleHemisphere samples a random direction from the canonical hemisphere whose normal is (0,0,1)
 * @return tuple of {sampledDir, probDensity} where sampledDir is a Vector3f in the canonical z-up space, and probDensity is the corresponding pdf val.
 */
std::tuple<Eigen::Vector3f, float> uniformSampleHemisphere() {
    float rand1 = generateRandFloat01();
    float rand2 = generateRandFloat01();
    float phi = 2 * M_PI * rand1;
    float theta = std::acos(1-rand2);
    float probDensity = 1.f / (2*M_PI); // == pdf(dir).

    return std::make_tuple(
                sphericalToCartesian(theta, phi),
                probDensity
            );
}

/**
 * @brief importanceSampleGlossyBRDF samples a direction from the canonical hemisphere using a PDF proportional to cos^n(angle between mirror direction and sample)
 * @param worldNormal normal at the intersection point in world space.
 * @param world_wi incident ray direction at the intersection point in world space. Must point away from the intersection point.
 * @param specularExp specular exponent (float) of the glossy surface defined in the scenefile
 * @return tuple of {sampledDir, probDensity} where sampledDir is a Vector3f in the canonical z-up space, and probDensity is the corresponding pdf val.
 */
std::tuple<Eigen::Vector3f, float> importanceSampleGlossyBRDF(Eigen::Vector3f &worldNormal, Eigen::Vector3f &world_wi, float specularExp) {
    float rand1 = generateRandFloat01();
    float rand2 = generateRandFloat01();
    Eigen::Vector3f z = Eigen::Vector3f::UnitZ();
    // glossy specular (phong) BRDF importance sampling: sample according to p = k*[cos(theta)]^n
    float phi = 2 * M_PI * rand1;
    float theta = std::acos(  pow(1 - rand2, (1.f/(specularExp+1)))  );
    Eigen::Vector3f sampledDir = sphericalToCartesian(theta, phi);
    float probDensity = ((specularExp+1) / (2*M_PI)) * pow(sampledDir.dot(z), specularExp); // p = k*[cos(theta)]^n
    // get the reflected/mirror direction in the canonical z-up space
    Eigen::Vector3f reflectedDirWorldSpace = reflectAboutNormal(worldNormal, world_wi);
    Eigen::Vector3f reflectedDirCanonical = rotationTo(worldNormal, z) * reflectedDirWorldSpace;

    // rotate sampled dir from canonical z-up space to to reflectedDirCanonical-up space
    sampledDir = rotationTo(z, reflectedDirCanonical) * sampledDir;

    return std::make_tuple(sampledDir, probDensity);
}

/**
 * @brief importanceSampleDiffuseBRDF samples a direction from the canonical hemisphere using a PDF proportional to cos(angle between normal and sample)
 * @return tuple of {sampledDir, probDensity} where sampledDir is a Vector3f in the canonical z-up space, and probDensity is the corresponding pdf val.
 */
std::tuple<Eigen::Vector3f, float> importanceSampleDiffuseBRDF() {
    float rand1 = generateRandFloat01();
    float rand2 = generateRandFloat01();
    // prob density is proportional to the angle between the normal and sampled dir. Scale factor comes from integrating p(w).
    float phi = 2 * M_PI * rand1;
    float theta = std::acos(1-rand2)/2.f; // the 1/2 scale factor comes from applying the inversion method with p=kcos(theta)
    Eigen::Vector3f sampledDirCanonical = sphericalToCartesian(theta, phi);
    float probDensity = (2.f/M_PI) * sampledDirCanonical.dot(Eigen::Vector3f::UnitZ()); // p = k*cos(theta) in the z-up space

    return std::make_tuple(sampledDirCanonical, probDensity);
}


/**
 * @param normal NORMALIZED world space normal at an intersection point
 * @param wi incident ray direction at the intersection point in world space. Must point away from the intersection point.
 * @param materialType the MaterialType of the intersected triangle
 * @param mat the Material of the intersected triangle
 * @return tuple of (1) a NORMALIZED direction vec sampled from a hemisphere whose base is perpendicular to the given normal,
 *              and (2) the probability density of the sample
 */
std::tuple<Eigen::Vector3f, float> sampleDirection(Eigen::Vector3f &normal, Eigen::Vector3f &wi, tinyobj::MaterialType materialType, const tinyobj::material_t &mat) {
    bool isDiffuse = (materialType == tinyobj::MaterialType::IDEAL_DIFFUSE);
    bool isGlossy = (materialType == tinyobj::MaterialType::GLOSSY_SPECULAR);
    // sample a direction and calculate its prob density
    Eigen::Vector3f sampledDirCanonical; // sample dir in z-up space, then transform to normal-up space
    float probDensity;
    Eigen::Vector3f z = Eigen::Vector3f::UnitZ();

    if (USE_BRDF_IMPORTANCE_SAMPLING && (isDiffuse || isGlossy)) {
        if (isDiffuse) {
            std::tie(sampledDirCanonical, probDensity) = importanceSampleDiffuseBRDF();
        } else {
            // glossy BRDF importance sampling
            std::tie(sampledDirCanonical, probDensity) = importanceSampleGlossyBRDF(normal, wi, mat.shininess);
        }
    } else {
        // uniform hemisphere sampling
        std::tie(sampledDirCanonical, probDensity) = uniformSampleHemisphere();
    }

    // rotate from canonical z-up space to normal-up space
    Eigen::Vector3f sampledDirWorld = rotationTo(z, normal) * sampledDirCanonical;

    return {sampledDirWorld.normalized(), probDensity};
}

Eigen::Vector2f sampleCirclePoint(float radius) {
    // uniformly sample polar coords
    float theta = 2 * M_PI * generateRandFloat01(); // in [0, 2pi]
    float r = radius * sqrt(generateRandFloat01()); // in [0, radius]
    // convert to Cartesian coords
    return Eigen::Vector2f(r*cos(theta), r*sin(theta));
}

float generateRandFloat01() {
    return uniformRealDist(generator);
}
