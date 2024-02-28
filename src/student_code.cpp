#include "student_code.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "halfEdgeMesh.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    vector<Vector2D> newPoints = std::vector<Vector2D>();
    for (int i = 0; i < points.size() - 1; i++) {
      Vector2D lerp = (1-t) * points[i] + (t) * points[i+1];
      newPoints.push_back(lerp);
    }
    return newPoints;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    vector<Vector3D> newPoints = std::vector<Vector3D>();
    for (int i = 0; i < points.size() - 1; i++) {
      Vector3D lerp = (1-t) * points[i] + (t) * points[i+1];
      newPoints.push_back(lerp);
    }
    return newPoints;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    vector<Vector3D> newPoints = evaluateStep(points, t);
    while (newPoints.size() > 1) {
      newPoints = evaluateStep(newPoints, t);
    }
    return newPoints[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {
    vector<Vector3D> rowPoints = std::vector<Vector3D>();  
    for (int i = 0; i < controlPoints.size(); i++) {
      rowPoints.push_back(evaluate1D(controlPoints[i], u));
    }
    Vector3D colPoints = evaluate1D(rowPoints, v);
    return colPoints;
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    HalfedgeCIter h = halfedge();      // get the outgoing half-edge of the vertex
    Vector3D normals = Vector3D();

    do {
        HalfedgeCIter h_twin = h->twin();
        HalfedgeCIter h_next = h_twin->next();
        VertexCIter v_twin = h_twin->vertex();
        VertexCIter v_next = h_next->twin()->vertex();
        Vector3D normal = cross(-v_twin->position, v_next->position);
        normals += normal;
        h = h_next;
    } while(h != halfedge());

    return normals / normals.norm();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    if (e0->isBoundary()) {
      return e0;
    }

    // Retrieve all elements
    HalfedgeIter bd = e0->halfedge();
    HalfedgeIter db = bd->twin();
    HalfedgeIter da = bd->next();
    HalfedgeIter ab = da->next();
    HalfedgeIter bc = db->next();
    HalfedgeIter cd = bc->next();

    VertexIter a = ab->vertex();
    VertexIter b = bd->vertex();
    VertexIter c = cd->vertex();
    VertexIter d = db->vertex();

    FaceIter abd = bd->face();
    FaceIter cbd = db->face();

    // Create new pointer names for the flip
    HalfedgeIter ac = bd;
    HalfedgeIter ca = db;

    // Set vertices
    a->halfedge() = ac;
    c->halfedge() = ca;

    // Set edges
    e0->halfedge() = ac;

    // Set faces
    FaceIter acd = abd;
    FaceIter cab = cbd;
    acd->halfedge() = ac;
    cab->halfedge() = ca;


    // Set halfedges
    ac->setNeighbors(cd, ca, a, e0, acd);
    cd->next() = da;
    da->next() = ac;

    ca->setNeighbors(ab, ac, c, ca->edge(), cbd);
    ab->next() = bc;
    bc->next() = ca;

    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {

    HalfedgeIter h = e0->halfedge();
    if (e0->isBoundary()) {
      return h->vertex();
    }

    VertexIter a = h->next()->next()->vertex();
    VertexIter b = h->vertex();
    VertexIter c = h->twin()->vertex();
    VertexIter d = h->twin()->next()->next()->vertex();

    HalfedgeIter ab = h->next()->next();
    HalfedgeIter bd = h->twin()->next();
    HalfedgeIter dc = h->twin()->next()->next();
    HalfedgeIter ca = h->next();

    // Create new vertex
    Vector3D endpoint1 = h->vertex()->position;
    Vector3D endpoint2 = h->twin()->vertex()->position;
    Vector3D midpoint = (endpoint1 + endpoint2) / 2;
    VertexIter m = newVertex();
    m->position = midpoint;

    // Create edges
    EdgeIter ea = newEdge();
    EdgeIter eb = e0;
    EdgeIter ec = newEdge();
    EdgeIter ed = newEdge();

    ea->halfedge() = h;
    ec->halfedge() = h;
    ed->halfedge() = h;

    // Create faces
    FaceIter fa = newFace();
    FaceIter fb = h->face();
    FaceIter fc = newFace();
    FaceIter fd = h->twin()->face();

    fa->halfedge() = h;
    fc->halfedge() = h;

    // Create halfedges
    HalfedgeIter ha = newHalfedge();
    HalfedgeIter hb = h;
    HalfedgeIter hc = newHalfedge();
    HalfedgeIter hd = newHalfedge();

    HalfedgeIter hat = newHalfedge();
    HalfedgeIter hbt = h->twin();
    HalfedgeIter hct = newHalfedge();
    HalfedgeIter hdt = newHalfedge();

    // Set halfedges
    ha->vertex() = a;
    hc->vertex() = c;
    hd->vertex() = d;

    hat->setNeighbors(ab, ha, m, ea, fb);
    hbt->vertex() = m;
    hbt->face() = fd;
    hct->setNeighbors(ca, hc, m, ec, fa);
    hdt->setNeighbors(dc, hd, m, ed, fc);

    ha->next() = hct;
    hb->next() = hat;
    hc->next() = hdt;
    hd->next() = hbt;

    ha->twin() = hat;
    hb->twin() = hbt;
    hc->twin() = hct;
    hd->twin() = hdt;

    ha->edge() = ea;
    hc->edge() = ec;
    hd->edge() = ed;

    ha->face() = fa;
    hc->face() = fc;
    hd->face() = fd;

    ca->next() = ha;
    ab->next() = hb;
    bd->next() = hd;
    dc->next() = dc;

    // Set faces
    fa->halfedge() = ha;
    // fc->halfedge() = hc;
    // fd->halfedge() = hbt;

    // Set edges
    ea->halfedge() = ha;
    ec->halfedge() = hc;
    ed->halfedge() = hd;
  
    return m;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.

  }
}
