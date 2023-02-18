#ifndef PATHTRACER_H
#define PATHTRACER_H

#include <QImage>

#include "scene/scene.h"


class PathTracer
{
public:
    PathTracer(int width, int height);

    void traceScene(QRgb *imageData, const Scene &scene);

private:
    int m_width, m_height;

    void toneMap(QRgb *imageData, std::vector<Eigen::Vector3f> &intensityValues);

    Eigen::Vector3f tracePixel(int x, int y, const Scene& scene, const Eigen::Matrix4f &invViewMatrix, int currentSampleIdx, int numSamples);

    Eigen::Vector3f traceRay(const Ray& r, const Scene &scene, bool countEmitted);

    Ray sampleRandomRayThruPixel(Eigen::Vector3f camOrigin, int imgRow, int imgCol);
    Ray stratifiedSubpixelSampling(Eigen::Vector3f camOrigin, int imgRow, int imgCol, int numCellsPerSide, int pixelCellIdx);

    float getContinuationProb();

    Eigen::Vector3f directLighting(IntersectionInfo& surfacePoint, const Eigen::Vector3f& w_i, const Scene& scene);

};

#endif // PATHTRACER_H
