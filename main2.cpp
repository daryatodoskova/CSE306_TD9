#include <iostream>
#include <vector>
#include <map>
#include <cmath>
// #define STB_IMAGE_IMPLEMENTATION
// #define STB_IMAGE_WRITE_IMPLEMENTATION
// #include "stb_image_write.h"
// #include "stb_image.h"


class Vector{
public:

    explicit Vector(double x = 0., double y = 0., double z = 0.) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;       
    }
    Vector operator+=(const Vector& b){
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
        return *this;
    }
    Vector operator*=(const Vector& b){
        data[0] *= b[0];
        data[1] *= b[1];
        data[2] *= b[2];
        return *this;
    }
    Vector operator-=(const Vector& b){
        data[0] -= b[0];
        data[1] -= b[1];
        data[2] -= b[2];
        return *this;
    }
    Vector operator/=(const Vector& b){
        data[0] /= b[0];
        data[1] /= b[1];
        data[2] /= b[2];
        return *this;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };

    
    int get_longest(){
        int res = 0;
        double maxi = data[0];
        for (int i = 1; i < 3; i++){
            if (data[i] > maxi){
                res = i;
                maxi = data[i];
            }
        }
        return res;
    }

private:
    double data[3];
};


Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const Vector& b) {
    return Vector(a[0] / b[0], a[1] / b[1], a[2] / b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}


Vector random_cos(const Vector& N){
    double r1 = ((double) rand() / RAND_MAX);
    double r2 = ((double) rand() / RAND_MAX);
    double x = sqrt(1 - r1) * cos(2. * M_PI * r2);
    double y = sqrt(1 - r1) * sin(2. * M_PI * r2);
    double z = sqrt(r1);

    Vector T1;
    if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
        T1 = Vector(0, N[2], -N[1]);
    } else {
        if (std::abs(N[1]) <= std::abs(N[0]) && std::abs(N[1]) <= std::abs(N[2])){
            T1 = Vector(N[2], 0, -N[0]);
        } else {
            T1 = Vector(N[1], -N[0], 0);
        }
    }
    T1.normalize();
    Vector T2 = cross(N, T1);

    return x*T1 + y*T2 + z*N;
}

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class TriangleMesh {
public:
    TriangleMesh() {};
    ~TriangleMesh() {};


    void center_scale(){

        Vector center(0,0,0);
        Vector minBbox(1E9, 1E9, 1E9);
        Vector maxBbox(-1E9, -1E9, -1E9);

        for (int i=0; i < vertices.size(); i++) {
            center += vertices[i];
            for (int j=0; j<3; j++) {
                minBbox[j] = std::min(minBbox[j], vertices[i][j]);
                maxBbox[j] = std::max(maxBbox[j], vertices[i][j]);
            }
        }

        center = center/vertices.size();
        double diagLength = (maxBbox - minBbox).norm();

        for (int i=0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] - center + Vector(0.5, 0.5, 0);
            vertices[i] = vertices[i] / diagLength;
        }
    }

   void writeOBJ(const char* obj) {

		FILE* f = fopen(obj, "w+");

		for (int i = 0; i < vertices.size(); i++) {
			fprintf(f, "v %f %f %f\n", vertices[i][0], vertices[i][1], vertices[i][2]);
		}

		for (int i = 0; i < indices.size(); i++) {
			fprintf(f, "f %u %u %u\n", indices[i].vtxi + 1, indices[i].vtxj + 1, indices[i].vtxk + 1);
		}

		fclose(f);

	}

    void scale_and_translate(double scaling, const Vector& translation);

    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
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
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
};

void TriangleMesh::scale_and_translate(double scaling, const Vector& translation) {
    for (int i=0; i < this->vertices.size(); i++) {
        this->vertices[i] = this->vertices[i] * scaling + translation;
    } 
}

class Edge{
public:
    Edge(int a=0, int b=0): a(a), b(b)
    int a, b;

};

bool operator<(const Edge& e1, const Edge& e2) {
    int mine1 = std::min(e1.a, e1.b); 
    int maxe1 = std::max(e1.a, e1.b);
    int mine2 = std::min(e2.a, e2.b);
    int maxe2 = std::max(e2.a, e2.b);
    return std::pair<int, int>(mine1, maxe1) <std::pair<int, int>(mine2,maxe2);
}

int main() {
    TriangleMesh mesh;
    mesh.readOBJ("goethe.obj");
    mesh.center_scale();
    std::map<Edge, std::vector<int> > edge_to_tri;
    std::map<int, std::vector<int> > vtx_to_tri;

    for (int i=0; i < mesh.indices.size(); i++) {
        int vtx1 = mesh.indices[i].vtxi;
        int vtx2 = mesh.indices[i].vtxj;
        int vtx3 = mesh.indices[i].vtxk;

        edge_to_tri[Edge(vtx1, vtx2)].push_back(i);
        edge_to_tri[Edge(vtx3, vtx1)].push_back(i);
        edge_to_tri[Edge(vtx2, vtx3)].push_back(i);

        vtx_to_tri[vtx1].push_back(i);
        vtx_to_tri[vtx2].push_back(i);
        vtx_to_tri[vtx3].push_back(i);
    }

    std::vector<bool> is_vtx_on_boundary(mesh.vertices.size(), false);

    std::vector<Edge> boundary_edges;
    for (auto it = edge_to_tri.begin(); it != edge_to_tri.end(); ++it) {
        if (it->second.size() == 1){
            boundary_edges.push_back(it->first);
            is_vtx_on_boundary[it->first.a]  = true; 
            is_vtx_on_boundary[it->first.b]  = true; 
        }
    }

    std::vector<Edge> ordered_boundary_edges(boundary_edges.size());
    int cur_edge = 0;
    ordered_boundary_edges[0] = boundary_edges[0];

    for (int i = 1; i < boundary_edges.size(); i++)
    {
        for (int j = 0; j < boundary_edges.size(); j++)
        {
            if (boundary_edges[j].a == ordered_boundary_edges[i-1].b)
            {
                ordered_boundary_edges[i] =  boundary_edges[j];
                break;
            }
        }
    }

    for (int i = 0; i < ordered_boundary_edges.size(); i++)
    {
        double theta = i/(double)ordered_boundary_edges.size() * 2. * M_PI;
        Vector circle_vtx;
        circle_vtx[0] = cos(theta)*0.5 + 0.5;
        circle_vtx[1] = sin(theta)*0.5 + 0.5;
        circle_vtx[2] = 0;
        mesh.vertices[ordered_boundary_edges[i].a] = circle_vtx;
    }
    
    for (int iter = 0; iter < 5000; iter++) 
    {
        std::vector<Vector> updated_vertices(mesh.vertices.size());

        updated_vertices = mesh.vertices;

        for (int i = 0; i < mesh.vertices.size(); i++)
        {
            if (is_vtx_on_boundary[i]) continue;
            
            Vector avg_neighbors(0,0,0);
            int nb_neighbors=0;

            for (int j = 0; j < vtx_to_tri[i].size(); j++)
            {
                if (mesh.indices[vtx_to_tri[i][j]].vtxi != i)
                avg_neighbors += mesh.vertices[mesh.indices[vtx_to_tri[i][j]].vtxi];
                if (mesh.indices[vtx_to_tri[i][j]].vtxj != i)
                avg_neighbors += mesh.vertices[mesh.indices[vtx_to_tri[i][j]].vtxj];
                if (mesh.indices[vtx_to_tri[i][j]].vtxk != i)
                avg_neighbors += mesh.vertices[mesh.indices[vtx_to_tri[i][j]].vtxk];
                nb_neighbors += 2;
            }
            updated_vertices[i] = avg_neighbors/nb_neighbors;
        }
        mesh.vertices = updated_vertices;
    }
    mesh.writeOBJ("flattened.obj");
}
