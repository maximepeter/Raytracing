#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <ctime>
#include <fstream>
#include <iostream>
#define M_PI 3.14159265358979323846

#include <random>

#include <string>
#include <stdio.h>
#include <list>

// lien vers canva3D https://www.cadnav.com/ telechargement impossible
// https://free3d.com/fr/ marche mieux

static std::default_random_engine engine(10);                // random seed = 10 // generateur de nombre aleatoire
static std::uniform_real_distribution<double> uniform(0, 1); // generation d'une loi uniforme entre 0 et 1

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0)
    {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }
    double operator[](int i) const { return coords[i]; };
    double &operator[](int i) { return coords[i]; };
    double sqrNorm()
    {
        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
    }
    Vector get_normalized()
    {
        double n = sqrt(sqrNorm());
        return Vector(coords[0] / n, coords[1] / n, coords[2] / n);
    }

    Vector &operator+=(const Vector &a)
    {
        coords[0] += a[0];
        coords[1] += a[1];
        coords[2] += a[2];
        return *this;
    }

private:
    double coords[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector &a)
{
    return Vector(-a[0], -a[1], -a[2]);
}

Vector operator*(double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector &a, double b)
{
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double sqr(double x)
{
    return x * x;
}

Vector random_cos(const Vector &N)
{
    double u1 = uniform(engine);
    double u2 = uniform(engine);
    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vector T1;
    if (N[0] < N[1] && N[0] < N[2])
    {
        T1 = Vector(0, N[2], -N[1]);
    }
    else
    {
        if (N[1] < N[2] && N[1] < N[0])
        {
            T1 = Vector(N[2], 0, -N[0]);
        }
        else
        {
            T1 = Vector(N[1], -N[0], 0);
        }
    }
    T1 = T1.get_normalized();
    Vector T2 = cross(N, T1);
    return z * N + x * T1 + y * T2;
}

class Ray
{
public:
    Ray(const Vector &C, const Vector &u) : C(C), u(u)
    {
    }

    Vector C, u;
};

class Object
{
public:
    Object(){};
    virtual bool intersect(const Ray &r, Vector &P, Vector &normale, double &t, Vector &color) = 0;

    Vector albedo;
    bool isMirror, isTransparent;
};

class BoudingBox
{
public:
    bool intersect(const Ray &r)
    {
        //intersection avec les plans verticaux
        double t1x = (mini[0] - r.C[0]) / r.u[0],
               t2x = (maxi[0] - r.C[0]) / r.u[0];
        double txMin = std::min(t1x, t2x), txMax = std::max(t1x, t2x);

        //intersection avec les plans horizontaux
        double t1y = (mini[1] - r.C[1]) / r.u[1],
               t2y = (maxi[1] - r.C[1]) / r.u[1];
        double tyMin = std::min(t1y, t2y), tyMax = std::max(t1y, t2y);

        //intersection avec les plans 3eme dimensions
        double t1z = (mini[2] - r.C[2]) / r.u[2],
               t2z = (maxi[2] - r.C[2]) / r.u[2];
        double tzMin = std::min(t1z, t2z), tzMax = std::max(t1z, t2z);

        //max et min des min et max
        double tMax = std::min(txMax, std::min(tyMax, tzMax)),
               tMin = std::max(txMin, std::max(tyMin, tzMin));
        if (tMax < 0)
            return false;
        return tMax > tMin;
    }
    Vector mini, maxi;
};

class Noeud
{
public:
    Noeud *fg, *fd;
    BoudingBox b;
    int debut, fin;
};

class TriangleIndices
{
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;    // indices within the uv coordinates array
    int ni, nj, nk;       // indices within the normals array
    int group;            // face group
};

class TriangleMesh : public Object
{
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vector &albedo, bool mirror = false, bool transp = false)
    {
        this->albedo = albedo;
        isMirror = mirror;
        isTransparent = transp;
        BVH = new Noeud;
    };

    BoudingBox buildBB(int debut, int fin) // debut et fin : indice de TRIANGLES
    {
        BoudingBox bb;
        bb.mini = Vector(1E9, 1E9, 1E9);
        bb.maxi = Vector(-1E9, -1E9, -1E9);
        // for (int i = 0; i < vertices.size(); i++)
        // {
        for (int i = debut; i < fin; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // bb.mini[j] = std::min(bb.mini[j], vertices[i][j]);
                // bb.maxi[j] = std::max(bb.maxi[j], vertices[i][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxi][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxi][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxj][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxj][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxk][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxk][j]);
            }
        }
        return bb;
    }

    void buildBVH(Noeud *n, int debut, int fin)
    {
        n->debut = debut;
        n->fin = fin;
        n->b = buildBB(n->debut, n->fin);

        Vector diag = n->b.maxi - n->b.mini;
        int dim; // dimension a diviser en 2
        if (diag[0] >= diag[1] && diag[0] >= diag[2])
        {
            dim = 0;
        }
        else
        {
            if (diag[1] >= diag[0] && diag[1] >= diag[2])
            {
                dim = 1;
            }
            else
            {
                dim = 2;
            }
        }

        double milieu = (n->b.mini[dim] + n->b.maxi[dim]) * 0.5; //milieu de la boite englobante selon la direction dim
        int indice_pivot = n->debut;
        for (int i = n->debut; i < n->fin; i++)
        {
            double milieu_triangle = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim]) / 3;
            if (milieu_triangle < milieu)
            {
                std::swap(indices[i], indices[indice_pivot]);
                indice_pivot++;
            }
        }

        n->fg = NULL;
        n->fd = NULL;

        if (indice_pivot == debut || indice_pivot == fin || fin - debut < 5)
            return;

        n->fg = new Noeud;
        n->fd = new Noeud;

        buildBVH(n->fg, n->debut, indice_pivot);
        buildBVH(n->fd, indice_pivot, n->fin);
    }

    void readOBJ(const char *obj)
    {

        char matfile[255];
        char grp[255];

        FILE *f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f))
        {
            char line[255];
            if (!fgets(line, 255, f))
                break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's')
            {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ')
            {
                Vector vec;

                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
                {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
                }
                else
                {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n')
            {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't')
            {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f')
            {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char *consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9)
                {
                    if (i0 < 0)
                        t.vtxi = vertices.size() + i0;
                    else
                        t.vtxi = i0 - 1;
                    if (i1 < 0)
                        t.vtxj = vertices.size() + i1;
                    else
                        t.vtxj = i1 - 1;
                    if (i2 < 0)
                        t.vtxk = vertices.size() + i2;
                    else
                        t.vtxk = i2 - 1;
                    if (j0 < 0)
                        t.uvi = uvs.size() + j0;
                    else
                        t.uvi = j0 - 1;
                    if (j1 < 0)
                        t.uvj = uvs.size() + j1;
                    else
                        t.uvj = j1 - 1;
                    if (j2 < 0)
                        t.uvk = uvs.size() + j2;
                    else
                        t.uvk = j2 - 1;
                    if (k0 < 0)
                        t.ni = normals.size() + k0;
                    else
                        t.ni = k0 - 1;
                    if (k1 < 0)
                        t.nj = normals.size() + k1;
                    else
                        t.nj = k1 - 1;
                    if (k2 < 0)
                        t.nk = normals.size() + k2;
                    else
                        t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else
                {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6)
                    {
                        if (i0 < 0)
                            t.vtxi = vertices.size() + i0;
                        else
                            t.vtxi = i0 - 1;
                        if (i1 < 0)
                            t.vtxj = vertices.size() + i1;
                        else
                            t.vtxj = i1 - 1;
                        if (i2 < 0)
                            t.vtxk = vertices.size() + i2;
                        else
                            t.vtxk = i2 - 1;
                        if (j0 < 0)
                            t.uvi = uvs.size() + j0;
                        else
                            t.uvi = j0 - 1;
                        if (j1 < 0)
                            t.uvj = uvs.size() + j1;
                        else
                            t.uvj = j1 - 1;
                        if (j2 < 0)
                            t.uvk = uvs.size() + j2;
                        else
                            t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3)
                        {
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0)
                                t.vtxi = vertices.size() + i0;
                            else
                                t.vtxi = i0 - 1;
                            if (i1 < 0)
                                t.vtxj = vertices.size() + i1;
                            else
                                t.vtxj = i1 - 1;
                            if (i2 < 0)
                                t.vtxk = vertices.size() + i2;
                            else
                                t.vtxk = i2 - 1;
                            if (k0 < 0)
                                t.ni = normals.size() + k0;
                            else
                                t.ni = k0 - 1;
                            if (k1 < 0)
                                t.nj = normals.size() + k1;
                            else
                                t.nj = k1 - 1;
                            if (k2 < 0)
                                t.nk = normals.size() + k2;
                            else
                                t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true)
                {
                    if (consumedline[0] == '\n')
                        break;
                    if (consumedline[0] == '\0')
                        break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3)
                    {
                        if (i0 < 0)
                            t2.vtxi = vertices.size() + i0;
                        else
                            t2.vtxi = i0 - 1;
                        if (i2 < 0)
                            t2.vtxj = vertices.size() + i2;
                        else
                            t2.vtxj = i2 - 1;
                        if (i3 < 0)
                            t2.vtxk = vertices.size() + i3;
                        else
                            t2.vtxk = i3 - 1;
                        if (j0 < 0)
                            t2.uvi = uvs.size() + j0;
                        else
                            t2.uvi = j0 - 1;
                        if (j2 < 0)
                            t2.uvj = uvs.size() + j2;
                        else
                            t2.uvj = j2 - 1;
                        if (j3 < 0)
                            t2.uvk = uvs.size() + j3;
                        else
                            t2.uvk = j3 - 1;
                        if (k0 < 0)
                            t2.ni = normals.size() + k0;
                        else
                            t2.ni = k0 - 1;
                        if (k2 < 0)
                            t2.nj = normals.size() + k2;
                        else
                            t2.nj = k2 - 1;
                        if (k3 < 0)
                            t2.nk = normals.size() + k3;
                        else
                            t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else
                    {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2)
                        {
                            if (i0 < 0)
                                t2.vtxi = vertices.size() + i0;
                            else
                                t2.vtxi = i0 - 1;
                            if (i2 < 0)
                                t2.vtxj = vertices.size() + i2;
                            else
                                t2.vtxj = i2 - 1;
                            if (i3 < 0)
                                t2.vtxk = vertices.size() + i3;
                            else
                                t2.vtxk = i3 - 1;
                            if (j0 < 0)
                                t2.uvi = uvs.size() + j0;
                            else
                                t2.uvi = j0 - 1;
                            if (j2 < 0)
                                t2.uvj = uvs.size() + j2;
                            else
                                t2.uvj = j2 - 1;
                            if (j3 < 0)
                                t2.uvk = uvs.size() + j3;
                            else
                                t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else
                        {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2)
                            {
                                if (i0 < 0)
                                    t2.vtxi = vertices.size() + i0;
                                else
                                    t2.vtxi = i0 - 1;
                                if (i2 < 0)
                                    t2.vtxj = vertices.size() + i2;
                                else
                                    t2.vtxj = i2 - 1;
                                if (i3 < 0)
                                    t2.vtxk = vertices.size() + i3;
                                else
                                    t2.vtxk = i3 - 1;
                                if (k0 < 0)
                                    t2.ni = normals.size() + k0;
                                else
                                    t2.ni = k0 - 1;
                                if (k2 < 0)
                                    t2.nj = normals.size() + k2;
                                else
                                    t2.nj = k2 - 1;
                                if (k3 < 0)
                                    t2.nk = normals.size() + k3;
                                else
                                    t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else
                            {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1)
                                {
                                    if (i0 < 0)
                                        t2.vtxi = vertices.size() + i0;
                                    else
                                        t2.vtxi = i0 - 1;
                                    if (i2 < 0)
                                        t2.vtxj = vertices.size() + i2;
                                    else
                                        t2.vtxj = i2 - 1;
                                    if (i3 < 0)
                                        t2.vtxk = vertices.size() + i3;
                                    else
                                        t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else
                                {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        fclose(f);
    }

    bool intersect(const Ray &r, Vector &P, Vector &normale, double &t, Vector &color)
    {
        if (!BVH->b.intersect(r))
            return false;

        t = 1E9;
        bool has_inter = false;

        std::list<Noeud *> l;
        l.push_back(BVH);
        while (!l.empty())
        {
            Noeud *c = l.front();
            l.pop_front();
            if (c->fg)
            {
                if (c->fg->b.intersect(r))
                {
                    l.push_front(c->fg);
                }
                if (c->fd->b.intersect(r))
                {
                    l.push_front(c->fd);
                }
            }
            else
            {

                for (int i = c->debut; i < c->fin; i++)
                {
                    // calcul d'intersection
                    const Vector &A = vertices[indices[i].vtxi],
                                 &B = vertices[indices[i].vtxj],
                                 &C = vertices[indices[i].vtxk];
                    Vector e1 = B - A,
                           e2 = C - A,
                           N = cross(e1, e2),
                           AO = r.C - A,
                           AOu = cross(AO, r.u);
                    double invUN = 1. / dot(r.u, N);

                    double beta = -dot(e2, AOu) * invUN,
                           gamma = dot(e1, AOu) * invUN,
                           alpha = 1 - beta - gamma,
                           localt = -dot(AO, N) * invUN;
                    if (beta >= 0 && gamma >= 0 && beta <= 1 && gamma <= 1 && alpha >= 0 && localt > 0)
                    {
                        has_inter = true;
                        if (localt < t)
                        {
                            t = localt;
                            // normale = N.get_normalized();
                            normale = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
                            normale = normale.get_normalized();
                            P = r.C + r.u;
                            int H = Htex[indices[i].group],
                                W = Wtex[indices[i].group];
                            Vector UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk];
                            UV = UV * Vector(W, H, 0);
                            int uvx = UV[0] + 0.5;
                            int uvy = UV[1] + 0.5;
                            uvx = uvx % W;
                            uvy = uvy % H;
                            if (uvx < 0)
                                uvx += W;
                            if (uvy < 0)
                                uvy += H;
                            uvy = H - uvy - 1;
                            color = Vector(std::pow(textures[indices[i].group][(uvy * W + uvx) * 3] / 255., 2.2),
                                           std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 1] / 255., 2.2),
                                           std::pow(textures[indices[i].group][(uvy * W + uvx) * 3 + 2] / 255., 2.2));
                        }
                    }
                }
            }
        }

        // for (int i = 0; i < indices.size(); i++)
        // {
        //     // calcul d'intersection
        //     const Vector &A = vertices[indices[i].vtxi],
        //                  &B = vertices[indices[i].vtxj],
        //                  &C = vertices[indices[i].vtxk];
        //     Vector e1 = B - A,
        //            e2 = C - A,
        //            N = cross(e1, e2),
        //            AO = r.C - A,
        //            AOu = cross(AO, r.u);
        //     double invUN = 1. / dot(r.u, N);

        //     double beta = -dot(e2, AOu) * invUN,
        //            gamma = dot(e1, AOu) * invUN,
        //            alpha = 1 - beta - gamma,
        //            localt = -dot(AO, N) * invUN;
        //     if (beta >= 0 && gamma >= 0 && beta <= 1 && gamma <= 1 && alpha >= 0 && localt > 0)
        //     {
        //         has_inter = true;
        //         if (localt < t)
        //         {
        //             t = localt;
        //             normale = N.get_normalized();
        //             P = r.C + r.u;
        //         }
        //     }
        // }

        return has_inter;
    };

    void loadTexture(const char *filename)
    {
        int W, H, C;
        unsigned char *texture = stbi_load(filename, &W, &H, &C, 3);
        Wtex.push_back(W);
        Htex.push_back(H);
        textures.push_back(texture);
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    std::vector<unsigned char *> textures;
    std::vector<int> Wtex, Htex;
    BoudingBox bb;

    Noeud *BVH;
};

class Sphere : public Object
{
public:
    Sphere(const Vector &O, double R, const Vector &albedo, bool isMirror = false, bool isTransparent = false) : O(O), R(R)
    {
        this->albedo = albedo;
        this->isMirror = isMirror;
        this->isTransparent = isTransparent;
    }
    bool intersect(const Ray &r, Vector &P, Vector &N, double &t, Vector &color)
    {
        // solves a*t^2 + b*t + c = 0
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C - O).sqrNorm() - R * R;
        double delta = b * b - 4 * a * c;
        color = this->albedo;

        if (delta < 0)
        {
            return false;
        }
        double sqDelta = sqrt(delta);
        double t2 = (-b + sqDelta) / (2 * a);

        if (t2 < 0)
            return false;

        // double t;
        double t1 = (-b - sqDelta) / (2 * a);

        if (t1 > 0)
        {
            t = t1;
        }
        else
        {
            t = t2;
        }

        P = r.C + t * r.u;
        N = (P - O).get_normalized();

        return true;
    }
    Vector O;
    double R;
};

class Scene
{
public:
    Scene(){};
    std::vector<Object *> objects;
    Vector L;
    double I;
    bool intersect(const Ray &r, Vector &P, Vector &N, Vector &albedo, bool &mirror, bool &transp, double &t, int &objectId)
    {
        t = 1E10;
        bool has_inter = false;
        for (int i = 0; i < objects.size(); i++)
        {
            Vector localP, localN, localAlbedo;
            double localt;
            if (objects[i]->intersect(r, localP, localN, localt, localAlbedo) && localt < t)
            {
                t = localt;
                has_inter = true;
                albedo = localAlbedo;
                P = localP;
                N = localN;
                mirror = objects[i]->isMirror;
                transp = objects[i]->isTransparent;
                objectId = i;
            }
        }
        return has_inter;
    };

    Vector getColor(const Ray &r, int rebond, bool lastDiffuse)
    {
        double epsilon = 0.00001;
        Vector P, N, albedo;
        double t;
        bool mirror, transp;
        int objectId;
        bool inter = intersect(r, P, N, albedo, mirror, transp, t, objectId);
        Vector color(0, 0, 0);
        if (rebond > 10) //5 normalement
            return Vector(0., 0., 0.);
        if (inter)
        {
            if (objectId == 0)
            {
                if (rebond == 0 || !lastDiffuse)
                {
                    return Vector(I, I, I) / (4 * M_PI * M_PI * sqr(dynamic_cast<Sphere *>(objects[0])->R));
                }
                return Vector(0., 0., 0.);
            }
            if (mirror)
            {
                Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                Ray reflectedRay(P + epsilon * N, reflectedDir);
                return getColor(reflectedRay, rebond + 1, false);
            }
            else
            {
                if (transp)
                {
                    double n1 = 1, n2 = 1.4;
                    Vector N2 = N;
                    if (dot(r.u, N) > 0)
                    { //on sort de la sphère
                        std::swap(n1, n2);
                        N2 = -N;
                    }

                    double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));

                    if (rad < 0) //rayon rasant, i-e rayon reflechie uniquement
                    {
                        Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                        Ray reflectedRay(P + epsilon * N, reflectedDir);
                        return getColor(reflectedRay, rebond + 1, false);
                    }

                    double k0 = sqr(n1 - n2) / sqr(n1 + n2);

                    double R = k0 + (1 - k0) * sqr(sqr(1 - std::abs(dot(N2, r.u)))) * (1 - std::abs(dot(N2, r.u)));

                    // implementation de la transmission de Fresnel avec tirage aleatoire et moyenne
                    double random = rand() % 100;

                    if (random <= 100 * R)
                    {
                        Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                        Ray reflectedRay(P + epsilon * N, reflectedDir);
                        return getColor(reflectedRay, rebond + 1, false);
                    }
                    else
                    {
                        // double T = 1 - R;
                        Vector tT = n1 / n2 * (r.u - dot(r.u, N2) * N2);
                        Vector tN = -sqrt(rad) * N2;
                        Vector refractedDir = tT + tN;
                        return getColor(Ray(P - epsilon * N2, refractedDir), rebond + 1, false);
                    }

                    // Vector reflectedDir = r.u - 2 * dot(r.u, N) * N;
                    // Ray reflectedRay(P + epsilon * N, reflectedDir);

                    // double T = 1 - R;

                    // Vector tT = n1 / n2 * (r.u - dot(r.u, N2) * N2);
                    // Vector tN = -sqrt(rad) * N2;
                    // Vector refractedDir = tT + tN;

                    // // sans transmission de Fresnel
                    // // return getColor(Ray(P - epsilon * N2, refractedDir), rebond + 1);

                    // // avec transmission de Fresnel
                    // return (T * getColor(Ray(P - epsilon * N2, refractedDir), rebond + 1) + R * getColor(reflectedRay, rebond + 1));
                }
                else //eclairage direct
                {
                    // Vector PL = L - P;
                    // double d = sqrt(PL.sqrNorm());
                    // Vector shadowP, shadowN, shadowAlbedo;
                    // double shadowt;
                    // int objectId;
                    // bool shadowMirror, shadowTransp;
                    // Ray shadowRay(P + 0.00001 * N, PL / d);
                    // bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt, objectId);
                    // if (shadowInter && shadowt < d)
                    // {
                    //     color = Vector(0., 0., 0.);
                    // }
                    // else
                    // {
                    //     color = (I / (4 * M_PI * d * d)) * (albedo / M_PI) * std::max(0., dot(N, PL / d));
                    // }

                    // eclairage direct

                    Vector PL = L - P;
                    PL = PL.get_normalized();
                    Vector w = random_cos(-PL);
                    Vector xprime = w * dynamic_cast<Sphere *>(objects[0])->R + dynamic_cast<Sphere *>(objects[0])->O;
                    Vector Pxprime = xprime - P;
                    double d = sqrt(Pxprime.sqrNorm());
                    Pxprime = Pxprime / d;

                    Vector shadowP, shadowN, shadowAlbedo;
                    double shadowt;
                    int objectId;
                    bool shadowMirror, shadowTransp;
                    Ray shadowRay(P + 0.00001 * N, Pxprime);
                    bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadowt, objectId);
                    if (shadowInter && shadowt < d - 0.0001)
                    {
                        color = Vector(0., 0., 0.);
                    }
                    else
                    {
                        double R2 = sqr(dynamic_cast<Sphere *>(objects[0])->R);
                        double proba = std::max(1E-8, dot(-PL, w)) / (M_PI * R2);
                        double J = std::max(0., dot(w, -Pxprime)) / (d * d);
                        color = (I / (4 * M_PI * M_PI * R2)) * (albedo / M_PI) * std::max(0., dot(N, Pxprime)) * J / proba;
                    }

                    // eclairage indirect
                    Vector wi = random_cos(N);
                    Ray wiRay(P + epsilon * N, wi);
                    color += albedo * getColor(wiRay, rebond + 1, true);
                }
            }
        }
        return color;
    }
};

void integrateCos()
{
    int N = 10000;
    double sigma = 0.25;
    double s = 0;
    for (int i = 0; i < N; i++)
    {
        double u1 = uniform(engine), u2 = uniform(engine);
        double xi = sigma * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
        if (xi >= -M_PI / 2 && xi <= M_PI / 2)
        {
            double p = 1 / (sigma * sqrt(2 * M_PI)) * exp(-xi * xi / (2 * sigma * sigma));
            s += pow(cos(xi), 10) / p / N;
        }
    }

    std::ofstream execution_file;
    execution_file.open("integral.txt");
    execution_file << s;
    execution_file.close();

    // std::cout << s << std::endl;
}

void integrate4D()
{
    int N = 100000;
    double sigma = 1;
    double s = 0;
    for (int i = 0; i < N; i++)
    {
        double u1 = uniform(engine), u2 = uniform(engine);
        double x1 = sigma * cos(2 * M_PI * u1) * sqrt(-2 * log(u2)),
               x2 = sigma * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));

        double u3 = uniform(engine), u4 = uniform(engine);
        double x3 = sigma * cos(2 * M_PI * u3) * sqrt(-2 * log(u4)),
               x4 = sigma * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

        if ((x1 >= -M_PI / 2 && x1 <= M_PI / 2) && (x2 >= -M_PI / 2 && x2 <= M_PI / 2) && (x3 >= -M_PI / 2 && x3 <= M_PI / 2) && (x4 >= -M_PI / 2 && x4 <= M_PI / 2))
        {
            double p1 = 1 / (sigma * sqrt(2 * M_PI)) * exp(-x1 * x1 / (2 * sigma * sigma)),
                   p2 = 1 / (sigma * sqrt(2 * M_PI)) * exp(-x2 * x2 / (2 * sigma * sigma)),
                   p3 = 1 / (sigma * sqrt(2 * M_PI)) * exp(-x3 * x3 / (2 * sigma * sigma)),
                   p4 = 1 / (sigma * sqrt(2 * M_PI)) * exp(-x4 * x4 / (2 * sigma * sigma));
            s += pow(cos(x1 + x2 + x3 + x4), 2) / (p1 * p2 * p3 * p4) / N;
        }
    }

    std::ofstream execution_file;
    execution_file.open("integral4D.txt");
    execution_file << s;
    execution_file.close();

    // std::cout << s << std::endl;
}

int main()
{
    float ini_time = clock();
    int W = 512;
    int H = 512;
    // integrateCos();
    // integrate4D();
    // return 0;

    Vector C(0, 0, 55);
    Scene scene;
    scene.I = 5E9;
    scene.L = Vector(-10, 70, 40);

    Sphere Slum(scene.L, 5, Vector(1., 1., 1.));
    Sphere S1(Vector(-17, 0, 0), 10, Vector(0., 0.5, 1.));
    Sphere S2(Vector(0, 10, 0), 10, Vector(1., 1., 1.), true);
    Sphere S3(Vector(17, 20, 0), 10, Vector(1., 0., 0.), false, true);
    Sphere Smurga(Vector(-1000, 0, 0), 970, Vector(0., 0., 1.));
    Sphere Smurdr(Vector(1000, 0, 0), 970, Vector(1., 0., 0.));
    Sphere Smurfa(Vector(0, 0, -1000), 940, Vector(0., 1., 0.));
    Sphere Smurde(Vector(0, 0, 1000), 940, Vector(1., 0., 1.));
    Sphere Ssol(Vector(0, -1000, 0), 990, Vector(1., 1., 1.), false);
    Sphere Splafond(Vector(0, 1000, 0), 990, Vector(1., 1., 1.));
    //TriangleMesh m(Vector(0., 1., 1.), false, false);
    //TriangleMesh m2(Vector(1., 1., 1.), false, false);
    //TriangleMesh mtable(Vector(1., 1., 1.), false, false);
    Sphere SMm(Vector(20, 20, -10), 10, Vector(1., 1., 1.), true);
    Sphere STm(Vector(0, 0, 10), 10, Vector(1., 1., 1.), false, true);
    //m.readOBJ("./chien/13463_Australian_Cattle_Dog_v3.obj");
    //m.loadTexture("./chien/Australian_Cattle_Dog_dif.jpg");

    //m2.readOBJ("./dumbell/10499_Dumbells_v1_L3.obj");
    //m2.loadTexture("./dumbell/10499_Dumbells_v1_diffuse.jpg");

    //mtable.readOBJ("./table/table.obj");
    //mtable.loadTexture("./table/render 1.jpg");
    // modifier les donnees pour rapeticer ou agrandir l'image
    //for (int i = 0; i < m2.vertices.size(); i++)
    {
        // O vers la droite 1 vers le haut 2 profondeur
        //m2.vertices[i][1] -= 5;
        // m2.vertices[i][1] += 5;
        //m2.vertices[i][2] = -m2.vertices[i][2];
        //m2.vertices[i][2] += 20;
    }
    // for (int i = 0; i < mtable.vertices.size(); i++)
    // {
    //     mtable.vertices[i][0] -= 25;
    //     mtable.vertices[i][2] = -mtable.vertices[i][2];
    // }
    //for (int i = 0; i < m.vertices.size(); i++)
    {

        // inversion y et z
        //std::swap(m.vertices[i][1], m.vertices[i][2]);
        // inversion x et z
        //std::swap(m.vertices[i][0], m.vertices[i][2]);
        //m.vertices[i][1] -= 10;
    }
    //for (int i = 0; i < m.vertices.size(); i++)
    {
        //std::swap(m.normals[i][1], m.normals[i][2]);
        //std::swap(m.normals[i][0], m.normals[i][2]);
        //m.normals[i] = -m.normals[i];
    }
    //m.buildBVH(m.BVH, 0, m.indices.size());
    //m2.buildBVH(m2.BVH, 0, m2.indices.size());
    // mtable.buildBVH(mtable.BVH, 0, mtable.indices.size());
    // m.buildBB();

    scene.objects.push_back(&Slum);
    scene.objects.push_back(&S1);
    scene.objects.push_back(&S2);
    scene.objects.push_back(&S3);
    scene.objects.push_back(&Smurga);
    scene.objects.push_back(&Smurdr);
    scene.objects.push_back(&Smurfa);
    // scene.objects.push_back(&Smurde);
    scene.objects.push_back(&Ssol);
    // scene.objects.push_back(&Splafond);
    //scene.objects.push_back(&m);
    //scene.objects.push_back(&m2);
    // scene.objects.push_back(&mtable);
    //scene.objects.push_back(&SMm);

    // scene.objects.push_back(&STm);

    double fov = 60 * M_PI / 180;

    int nbrays = 40;
    double angleVertical = 0 * M_PI / 180, angleHorizontal = 0 * M_PI / 180;
    Vector up(0, cos(angleVertical), sin(angleVertical));
    Vector right(cos(angleHorizontal), 0, sin(angleHorizontal));

    double up0 = up[0];
    up[0] = cos(angleHorizontal) * up[0] - sin(angleHorizontal) * up[2];
    up[2] = sin(angleHorizontal) * up0 + cos(angleHorizontal) * up[2];

    Vector viewDirection = cross(up, right);

    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {

            //2
            Vector color(0, 0, 0);
            for (int k = 0; k < nbrays; k++)
            {

                double u1 = uniform(engine), u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2)),
                       x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));

                double u3 = uniform(engine), u4 = uniform(engine);
                double x3 = 0.01 * cos(2 * M_PI * u3) * sqrt(-2 * log(u4)), // remettre à 1
                    x4 = 0.01 * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

                Vector u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, W / (2. * tan(fov / 2)));
                u = u.get_normalized();
                u = u[0] * right + u[1] * up + u[2] * viewDirection;

                Vector target = C + 55 * u;
                Vector Cprime = C + Vector(x3, x4, 0);
                Vector uprime = (target - Cprime).get_normalized();
                // Ray r(C, u);
                Ray r(Cprime, uprime);

                color += scene.getColor(r, 0, false);
            }

            color = color / nbrays;

            // implementation de la transmission de Fresnel avec tirage aleatoire et moyenne
            // Vector color = Vector(0., 0., 0.);
            // int nbVal = 30;
            // for (int i = 0; i < nbVal; i++)
            // {
            //     color = color + scene.getColor(r, 0);
            // }
            // color = color / nbVal;

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
        }
    }

    // wwriting the image on a png file
    stbi_write_png("raytracing.png", W, H, 3, &image[0], 0);

    // writing execution time on a txt file
    std::ofstream execution_file;
    execution_file.open("time.txt");
    execution_file << (clock() - ini_time) / CLOCKS_PER_SEC;
    execution_file.close();

    return 0;
};

// run command line
// g++ raytracing.cpp -o raytracing.out
// ./raytracing.out 