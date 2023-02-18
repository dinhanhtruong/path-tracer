#include "pathtracer.h"

#include <iostream>

#include <Eigen/Dense>

#include <util/CS123Common.h>

using namespace Eigen;
//std::default_random_engine generator;
//std::uniform_real_distribution<float> uniformDist(0.0,1.0);

PathTracer::PathTracer(int width, int height)
    : m_width(width), m_height(height)
{
}

void PathTracer::traceScene(QRgb *imageData, const Scene& scene) {
    std::vector<Vector3f> intensityValues(m_width * m_height);
    Matrix4f invViewMat = (scene.getCamera().getScaleMatrix() * scene.getCamera().getViewMatrix()).inverse();  // scale???
    int numSamplesPerPixel = USE_STRATIFIED_SUBPIXEL_SAMPLING ? pow(floor(sqrt(NUM_SAMPLES_PER_PIXEL)), 2) // get closest smaller perfect square
                                                              : NUM_SAMPLES_PER_PIXEL;
    #pragma omp parallel for
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            // trace N ray samples through each pixel and accumulate the radiance
            Vector3f currPixelRadiance = Vector3f::Zero();
            for(int i=0; i < numSamplesPerPixel; i++) {
                currPixelRadiance += tracePixel(x, y, scene, invViewMat, i, numSamplesPerPixel);
            }
            // store the averaged radiance
            int offset = x + (y * m_width);
            intensityValues[offset] = currPixelRadiance / numSamplesPerPixel;
            if (x == 0 && y % 20 == 0) {
                std::cout << (int)(100*y/m_height) << "%" << std::endl;
            }
        }
    }

    toneMap(imageData, intensityValues);
}

Ray PathTracer::sampleRandomRayThruPixel(Vector3f camOrigin, int imgRow, int imgCol) {
    // uniformly sample random point on the current pixel of the img plane by translating from the top left corner
    float randXOffset = generateRandFloat01(); // in [0,1]
    float randYOffset = -generateRandFloat01(); // in [-1,0]
    // direction is offset randomly from top left corner of pixel
    float imgPlaneDepth = DEPTH_OF_FIELD <= 0 ? 1 : DEPTH_OF_FIELD;
    Vector3f imagePlaneLoc(
                (2.f * (imgCol + randXOffset) / m_width) - 1,
                1 - (2.f * (imgRow + randYOffset) / m_height),
                -1 // must be on image plane (z=-1)
            );
    // rescale image plane based on its depth along look
    imagePlaneLoc *= imgPlaneDepth;
    Vector3f d = imagePlaneLoc - camOrigin;
    // construct camera-space ray
    return Ray(camOrigin, d.normalized());
}

Ray PathTracer::stratifiedSubpixelSampling(Vector3f camOrigin, int imgRow, int imgCol, int numCellsPerSide, int pixelCellIdx) {
    // basic idea: split pixel into numCellsPerSide*numCellsPerSide cells and sample one ray from each
    int cellRow = pixelCellIdx / numCellsPerSide;
    int cellCol = pixelCellIdx % numCellsPerSide;
    // assume pixels have width 1 => pixel cells have width 1/numCellsPerSide
    float cellWidth = 1.f/numCellsPerSide;
    float cellXOffset = cellWidth * cellCol;
    float cellYOffset = -cellWidth * cellRow;
    // randomly offset ray from the top left of the cell
    float randXOffset = generateRandFloat01() * cellWidth; // in [0,cellWidth]
    float randYOffset = -generateRandFloat01() * cellWidth; // in [-cellWidth,0]
    float imgPlaneDepth = DEPTH_OF_FIELD <= 0 ? 1 : DEPTH_OF_FIELD;
    Vector3f imagePlaneLoc(
                (2.f * (imgCol + cellXOffset + randXOffset)/m_width ) - 1,
                 1 - (2.f * (imgRow + cellYOffset + randYOffset)/m_height),
                -1 // must be on image plane (z=-1)
            );
    // rescale image plane based on its depth along look
    imagePlaneLoc *= imgPlaneDepth;
    Vector3f d = imagePlaneLoc - camOrigin;
    // construct camera-space ray
    return Ray(camOrigin, d.normalized());
}

Vector3f PathTracer::tracePixel(int x, int y, const Scene& scene, const Matrix4f &invViewMatrix, int currentSampleIdx, int numSamples) {
    Vector3f origin(0, 0, 0);
    if (DEPTH_OF_FIELD > 0) {
        // scatter the eye location over a circular planar lens normal to look
        Vector2f circlePoint = sampleCirclePoint(0.35);
        origin = Vector3f(circlePoint.x(), circlePoint.y(), 0);
    }

    Ray r(origin, Vector3f::Zero());
    // get camera space ray from the integer image plane coords
    if (USE_STRATIFIED_SUBPIXEL_SAMPLING) {
        // shoot one ray through each sub-pixel 'cell'
        r = stratifiedSubpixelSampling(origin, y, x, sqrt(numSamples), currentSampleIdx);
    } else {
        r = sampleRandomRayThruPixel(origin, y, x);
    }

    // transform ray to world space
    r = r.transform(invViewMatrix);
    return traceRay(r, scene, true); // always count direct lighting toward camera
}

/**
 * @brief PathTracer::traceRay
 * @param r world space ray
 * @param scene to be rendered
 * @return incoming radiance Li(r_o, r_d) toward the origin of the given ray r from its first intersection point
 */
Vector3f PathTracer::traceRay(const Ray& r, const Scene& scene, bool countEmitted) {
    IntersectionInfo insct;
    Ray rayOut(r);

    // find first intersection (if any) of ray w/ scene geometry
    if(scene.getIntersection(rayOut, &insct)) {

        Vector3f radianceIn(0,0,0); // incoming radiance Li(r_o, r_d) toward the origin of the ray r from the intersection

        const Triangle *tri = static_cast<const Triangle *>(insct.data); // cast intersected obj to tri since handling tri meshes
        Vector3f normal = tri->getNormal(insct.hit); // world space normal at insct point
        Vector3f wi_intersection = -rayOut.d; // incident ray direction at the insct point. Points outward from insct.

        const tinyobj::material_t mat = tri->getMaterial(); // get the material of the intersected mesh tri
        bool intersectionIsMirror = (tri->getMaterialType() == tinyobj::MaterialType::IDEAL_SPECULAR);
        bool intersectionIsRefractive = (tri->getMaterialType() == tinyobj::MaterialType::DIELECTRIC_REFRACTION);

        // if hit a light, accumulate emitted radiance from the light toward the ray origin
        if (tri->getMaterialType() == tinyobj::MaterialType::LIGHT_SOURCE) {
            if (countEmitted) {
                radianceIn += tri->emittedRadiance(wi_intersection);
                assert(radianceIn[0] > 0 && radianceIn[1] > 0 && radianceIn[2] > 0);
            }
            return radianceIn; // terminate path at light
        }
        if (!intersectionIsMirror && !intersectionIsRefractive){
            // count direct lighting to the intersection
            radianceIn += directLighting(insct, wi_intersection, scene);
        }

        // use single-sample MC estimate of the recursive integral term
        Vector3f nextDir; // world space
        float nextDirProb; // probability density of sampling the next ray direction
        float cosTheta; // theta = angle between normal and w_i at insct. Only needed for indirect lighting.
        if (intersectionIsMirror) {
            nextDir = reflectAboutNormal(normal, wi_intersection);
        } else if(intersectionIsRefractive) {
            // assume normal from getIntersection always faces outward
            bool enteringFromAir = (normal.dot(wi_intersection) > 0);
            // re-orient the normal so that it is in the same hemisphere as in the same hemisphere as wi_intersection
            normal = enteringFromAir ? normal : -normal;
            // determine index of refraction of the
            float iorIn;
            float iorOut;
            if (enteringFromAir) {
                // assume air has ior=1
                iorIn = 1;
                iorOut = mat.ior;
            } else {
                iorIn = mat.ior;
                iorOut = 1;
            }
            // get the refracted/transmitted ray (will be empty iff have total internal reflection)
            std::optional<Vector3f> refractedRay = refract(normal, wi_intersection, iorIn, iorOut);
            // special case: always reflect if have total internal reflection
            if (!refractedRay.has_value()) {
                nextDir = reflectAboutNormal(normal, wi_intersection);
            } else {
                // general case: use Schlick's approximation of the Fresnel term to determine the probability of reflection vs transmission/refraction.
                //               & use Snell's law to get the refracted ray dir
                float probReflect = schlickApprox(normal.dot(wi_intersection), iorIn, iorOut); // = prob(reflect)
                nextDir = (generateRandFloat01() < probReflect) ? reflectAboutNormal(normal, wi_intersection)
                                                                : refractedRay.value();
            }
        } else {
            // general case: sample from hemisphere for indirect lighting
            std::tie(nextDir, nextDirProb) = sampleDirection(normal, wi_intersection, tri->getMaterialType(), mat);
            cosTheta = normal.dot(nextDir);
        }

        Vector3f f_r = tri->brdf(nextDir, wi_intersection, normal); // material-dependent BRDF at intersection point
        float prob_rr = getContinuationProb(); // russian roulette
        if (generateRandFloat01() < prob_rr){
            // shoot recursive ray in the chosen dir
            Ray nextRay(insct.hit, nextDir);
            // compute incoming radiance integral. 2 cases: general indirect lighting & reflective/refractive surfaces
            if (intersectionIsMirror || intersectionIsRefractive) {
                // all incoming radiance is from the refracted/reflected direction: cos(theta) term cancels due to ideal BSDF.
                // approximate prob density = 1 and dirac delta = 1.
                radianceIn += Vector3f(
                                f_r.array() * traceRay(nextRay, scene, true).array()
                                /(prob_rr)
                            );
            } else {
                // add radiance contribution of the sampled ray
                // don't double count direct lighting (i.e. split integral) so that indirect + direct rays partition the hemisphere
                radianceIn += Vector3f(
                                f_r.array() * traceRay(nextRay, scene, false).array() * cosTheta
                                /(nextDirProb * prob_rr)
                            );
            }
        }
        return radianceIn;
    } else {
        return Vector3f(0, 0, 0);
    }
}

float PathTracer::getContinuationProb() {
    return 1.f - PATH_TERMINATION_PROB;
}

/**
 * @brief PathTracer::directLighting Samples a point on an area light source to determine the direct radiance from the light to p.
 *          0 if the sampled light point is occluded.
 * @param surfacePoint intersection point from which we want to calculate direct lighting
 * @param w_i incoming ray direction at p (outward facing)
 * @return radiance from the sampled light point p' to the surface point p
 */
Vector3f PathTracer::directLighting(IntersectionInfo& surfacePoint, const Vector3f& w_i, const Scene& scene) {
    Vector3f totalDirectRadiance = Vector3f::Zero();
    const std::vector<Triangle *>& lights = scene.getEmissives();
    Vector3f p(surfacePoint.hit);
    for (int i = 0; i < NUM_DIRECT_LIGHTING_SAMPLES; i++) {
        // sample rand emissive triangle
        int lightID = (int)(generateRandFloat01() * lights.size()) % lights.size();
        Triangle *light = lights[lightID];
        // sample point p' on triangle
        Vector3f lightPoint = light->sampleTrianglePoint();
        // trace ray from p to p'
        Vector3f vecToLight = lightPoint - p;
        IntersectionInfo lightInsct;
        Ray shadowRay(p, vecToLight.normalized());
        scene.getIntersection(shadowRay, &lightInsct); // shadowRay must intersect with the light XOR an occluding object
        // determine visibility of p: light is occluded if |traced ray| < |p - p'|   (i.e. intersection happened earlier)
        if (lightInsct.t < vecToLight.norm() - FLOAT_EPSILON) {
            continue; // no contribution to MC estimator
        }
        const Triangle *surface = static_cast<const Triangle *>(surfacePoint.data);
        Vector3f surfaceNormal = surface->getNormal(p); // world space normal at insct point
        Vector3f f_r = surface->brdf(shadowRay.d, w_i, surfaceNormal); // material-dependent BRDF
        Vector3f lightNormal = light->getNormal(lightPoint); // world space normal at insct point
        float lightArea = light->getTriangleArea();

        float cosTheta = (shadowRay.d).dot(surfaceNormal);
        float cosThetaLight = (-shadowRay.d).dot(lightNormal);
        // ignore contribution if angle between (surface/light) normal and shadow ray > 90 deg
        if (cosTheta < 0 || cosThetaLight < 0) continue;

        // add to MC estimator of radiance from light
        Vector3f currRadiance = lightArea * Vector3f( // multiply by light area <=> div by pdf=1/A
                                    f_r.array() *                                   // BRDF at the non-light surface from which we trace shadowRay
                                    light->emittedRadiance(-shadowRay.d).array() *  // L_o(p', w')
                                    (shadowRay.d).dot(surfaceNormal) *              // cos(theta)
                                    (-shadowRay.d).dot(lightNormal) /               // cos(theta')
                                    vecToLight.norm()                               // |p-p'|
                                );
        totalDirectRadiance += currRadiance;
    }
    return totalDirectRadiance / NUM_DIRECT_LIGHTING_SAMPLES;
}

void PathTracer::toneMap(QRgb *imageData, std::vector<Vector3f> &intensityValues) {
    // find the maximum luminance value in the entire image
    float L_max = 0;
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            Vector3f radiance = intensityValues[offset].cwiseMax(0);
            float curr_L = RGBradianceToLuminance(radiance);
            L_max = std::max(L_max, curr_L);
        }
    }
    // write to the output image
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);

            // clip to [0,infty]
            Vector3f radiance = intensityValues[offset].cwiseMax(0);
            if (radiance.isZero()) {
                imageData[offset] = qRgb(0, 0, 0); // avoid divide by 0 issues later
                continue;
            }

            // convert to luminance using convex combination of color channels
            float L_world = RGBradianceToLuminance(radiance);
            // apply (extended) reinhard tone map to luminance
            float numerator = L_world * (1.f + (L_world / (L_max * L_max)));
            float L_display = numerator / (1.f + L_world);
            // scale original radiance by luminance ratios
            Vector3f tonemappedRadiance = (L_display / L_world) * radiance;

            // gamma correction
            tonemappedRadiance = tonemappedRadiance.array().pow(1.f/2);
            // clip to [0,1]
            tonemappedRadiance = tonemappedRadiance.cwiseMin(1);

            // scale to [0,255] for displaying
            tonemappedRadiance *= 255.f;

            // cast to int RGB vals
            imageData[offset] = qRgb((int)tonemappedRadiance[0], (int)tonemappedRadiance[1], (int)tonemappedRadiance[2]);
        }
    }
}

