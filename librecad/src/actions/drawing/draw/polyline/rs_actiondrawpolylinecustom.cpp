#include "rs_actiondrawpolylinecustom.h"
#include "rs_actioninterface.h"
#include "rs_graphicview.h"
#include "rs_entitycontainer.h"
#include "rs_stringutils.h"
#include "rs_polyline.h"
#include "rs_line.h"
#include "rs_arc.h"
#include "rs_debug.h"

#include <QAction>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QSettings>
#include <cmath>

static const char* SETTINGS_GROUP = "RS_ActionDrawPolylineCustom";
static const char* SETTINGS_KEY_LAST_RADIUS = "lastRadius";

RS_ActionDrawPolylineCustom::RS_ActionDrawPolylineCustom(RS_EntityContainer* container,
                                                         RS_GraphicView* graphicView)
    : RS_ActionInterface("Draw polyline (CUSTOM)", container, graphicView),
      polyline(nullptr),
      status(SetStartPoint),
      radius(5.0), // default, will be overridden by settings
      variableRadius(false),
      settings()
{
    // load last radius
    settings.beginGroup(SETTINGS_GROUP);
    radius = settings.value(SETTINGS_KEY_LAST_RADIUS, radius).toDouble();
    settings.endGroup();
}

RS_ActionDrawPolylineCustom::~RS_ActionDrawPolylineCustom() {
    // store last used radius
    settings.beginGroup(SETTINGS_GROUP);
    settings.setValue(SETTINGS_KEY_LAST_RADIUS, radius);
    settings.endGroup();
}

void RS_ActionDrawPolylineCustom::init(int s) {
    RS_ActionInterface::init(s);
    status = (s==0) ? SetStartPoint : SetNextPoint;
    polyline = nullptr;
    graphicView->setMouseCursor(RS2::CadCursor);
    // display initial command prompt
    RS2::updateCommandWidget(tr("Specify first point"));
}

QAction* RS_ActionDrawPolylineCustom::createGUIAction(RS2::ActionType, QObject* parent) {
    QAction* action = new QAction(tr("Polyline (CUSTOM)"), parent);
    action->setToolTip(tr("Draw polyline and automatically fillet vertices with default radius"));
    return action;
}

void RS_ActionDrawPolylineCustom::trigger() {
    // Not used; geometry is created on clicks
}

void RS_ActionDrawPolylineCustom::mousePressEvent(QMouseEvent* e) {
    RS_DEBUG->print("RS_ActionDrawPolylineCustom::mousePressEvent");
    RS_Vector pos = snapPoint(e);
    if (e->button() == Qt::LeftButton) {
        if (status == SetStartPoint) {
            startPoint = pos;
            lastPoint = pos;
            // create polyline container entity (no vertices yet)
            polyline = new RS_Polyline(container, nullptr);
            // add first vertex
            polyline->addVertex(startPoint, 0.0); // bulge=0 -> straight
            container->addEntity(polyline);
            status = SetNextPoint;
            RS2::updateCommandWidget(tr("Specify next point or [Undo/Close]: R=%1").arg(radius));
        } else if (status == SetNextPoint) {
            // Add new segment from lastPoint to pos, and compute fillet at lastPoint (between prev and this)
            // If polyline has only one segment so far (only first vertex exists), we just add a straight vertex
            RS_Vector prev = lastPoint;
            RS_Vector newPt = pos;

            // If this is second segment (i.e., we already have at least 2 vertices),
            // compute fillet between prev segment and new segment and replace last vertex with tangent point + bulge
            int vcount = polyline->countVertices();
            if (vcount >= 2) {
                // get previous vertex before lastPoint
                RS_Vector prevPrev = polyline->getVertexAt(vcount-2);
                RS_Vector tangent1, tangent2, centerArc;
                double arcAngle;
                bool arcReversed;
                bool ok = computeFilletTangents(prevPrev, prev, newPt, radius, tangent1, tangent2, centerArc, arcAngle, arcReversed);

                if (ok) {
                    // replace the last vertex position with tangent1 (so the previous segment ends at tangent1)
                    polyline->setVertexAt(vcount-1, tangent1, 0.0);

                    // append tangent2 as the next straight vertex after the arc (bulge will encode arc between tangent1 and tangent2)
                    // compute bulge from arcAngle and direction
                    double includedAngle = fabs(arcAngle);
                    double bulge = tan(includedAngle / 4.0);
                    if (arcReversed) {
                        bulge = -bulge;
                    }
                    polyline->addVertex(tangent2, bulge);

                    // The final straight piece from tangent2 to newPt will be added next by just adding newPt as straight vertex.
                    polyline->addVertex(newPt, 0.0);

                    lastPoint = newPt;
                    RS2::updateCommandWidget(tr("Added vertex with fillet R=%1").arg(radius));
                    return;
                } else {
                    // cannot create fillet (too short or self-intersect), fallback to straight join
                    polyline->addVertex(newPt, 0.0);
                    lastPoint = newPt;
                    RS2::updateCommandWidget(tr("Added vertex (fillet skipped) R=%1").arg(radius));
                    return;
                }
            } else {
                // first segment only: add straight vertex
                polyline->addVertex(newPt, 0.0);
                lastPoint = newPt;
                RS2::updateCommandWidget(tr("Added vertex"));
                return;
            }
        }
    } else if (e->button() == Qt::RightButton) {
        // finish polyline
        if (status == SetNextPoint) {
            finishPolyline();
        } else {
            // cancel
            if (polyline) {
                container->deleteEntity(polyline);
                polyline = nullptr;
            }
            actionFinished();
        }
    }
}

void RS_ActionDrawPolylineCustom::mouseMoveEvent(QMouseEvent* e) {
    RS_Vector mouse = snapPoint(e);
    // dynamic preview: show line from lastPoint to mouse and fillet preview between previous and that
    graphicView->clearPreview();
    if (status == SetNextPoint && polyline) {
        // preview last straight segment
        RS_Line* previewLine = new RS_Line(nullptr, nullptr, lastPoint, mouse);
        graphicView->drawEntity(previewLine, RS2::FlagPreview);
        delete previewLine;

        // if we already have at least 2 vertices, preview fillet at lastPoint
        int vcount = polyline->countVertices();
        if (vcount >= 2) {
            RS_Vector prevPrev = polyline->getVertexAt(vcount-2);
            RS_Vector tangent1, tangent2, centerArc;
            double arcAngle;
            bool arcReversed;
            if (computeFilletTangents(prevPrev, lastPoint, mouse, radius, tangent1, tangent2, centerArc, arcAngle, arcReversed)) {
                // create arc preview
                double ang1 = (tangent1 - centerArc).angle();
                double ang2 = (tangent2 - centerArc).angle();
                RS_Arc* previewArc = new RS_Arc(nullptr, nullptr, centerArc, (tangent1-centerArc).magnitude(), ang1, ang2, arcReversed);
                graphicView->drawEntity(previewArc, RS2::FlagPreview);
                delete previewArc;

                // preview tangent lines (from prevPrev->tangent1 and tangent2->mouse)
                RS_Line* tline1 = new RS_Line(nullptr, nullptr, prevPrev, tangent1);
                RS_Line* tline2 = new RS_Line(nullptr, nullptr, tangent2, mouse);
                graphicView->drawEntity(tline1, RS2::FlagPreview);
                graphicView->drawEntity(tline2, RS2::FlagPreview);
                delete tline1;
                delete tline2;
            }
        }
    }
    graphicView->redraw();
}

void RS_ActionDrawPolylineCustom::mouseReleaseEvent(QMouseEvent* e) {
    RS_UNUSED(e);
}

void RS_ActionDrawPolylineCustom::keyPressEvent(QKeyEvent* e) {
    // Allow numeric entry + Enter to set radius (simple implementation)
    if (e->key() == Qt::Key_R) {
        // start typing radius or toggle variable mode
        variableRadius = !variableRadius;
        RS2::updateCommandWidget(tr("Variable radius %1").arg(variableRadius ? tr("ON") : tr("OFF")));
        return;
    }

    if (e->key() == Qt::Key_Escape) {
        // cancel
        if (polyline) {
            container->deleteEntity(polyline);
            polyline = nullptr;
        }
        actionFinished();
        return;
    }

    // accept digits, decimal point and backspace for radius input
    if ((e->key() >= Qt::Key_0 && e->key() <= Qt::Key_9) || e->key() == Qt::Key_Period || e->key() == Qt::Key_Backspace || e->key() == Qt::Key_Minus) {
        if (e->key() == Qt::Key_Backspace && !radiusTextBuffer.isEmpty()) {
            radiusTextBuffer.chop(1);
        } else if (e->text().length() > 0) {
            radiusTextBuffer += e->text();
        }
        RS2::updateCommandWidget(tr("Radius: %1").arg(radiusTextBuffer));
        return;
    }

    if (e->key() == Qt::Key_Return || e->key() == Qt::Key_Enter) {
        if (!radiusTextBuffer.isEmpty()) {
            bool ok=false;
            double v = RS_String::cleanNumber(radiusTextBuffer).toDouble(&ok);
            if (ok && v>=0.0) {
                setRadius(v);
            }
            radiusTextBuffer.clear();
        } else {
            // finish polyline
            if (status == SetNextPoint) finishPolyline();
        }
    }
}

void RS_ActionDrawPolylineCustom::commandEvent(RS_CommandEvent* e) {
    QString cmd = e->getCommand();
    // allow direct numeric radius: command like "r 10" or "10"
    bool ok=false;
    double val = RS_String::cleanNumber(cmd).toDouble(&ok);
    if (ok) {
        setRadius(val);
        RS2::updateCommandWidget(tr("Radius set to %1").arg(radius));
        return;
    }

    if (cmd.toLower() == "v" || cmd.toLower() == "var") {
        setVariableRadiusMode(!variableRadius);
        RS2::updateCommandWidget(tr("Variable radius %1").arg(variableRadius ? tr("ON") : tr("OFF")));
        return;
    }

    RS_ActionInterface::commandEvent(e);
}

void RS_ActionDrawPolylineCustom::mouseWheelEvent(QWheelEvent* e) {
    // Use Ctrl + wheel to change radius quickly, otherwise ignore
    if (e->modifiers() & Qt::ControlModifier) {
        int delta = e->angleDelta().y();
        double step = (e->modifiers() & Qt::ShiftModifier) ? 0.1 : 1.0;
        if (delta > 0) {
            setRadius(std::max(0.0, radius + step));
        } else {
            setRadius(std::max(0.0, radius - step));
        }
        RS2::updateCommandWidget(tr("Radius: %1").arg(radius));
    }
}

bool RS_ActionDrawPolylineCustom::updateMouseCursor() {
    graphicView->setMouseCursor(RS2::CadCursor);
    return true;
}

void RS_ActionDrawPolylineCustom::finishPolyline() {
    // finalize polyline: remove redundant consecutive identical vertices (if needed)
    if (!polyline) {
        actionFinished();
        return;
    }
    if (polyline->countVertices() < 2) {
        container->deleteEntity(polyline);
    } else {
        // optionally close if last vertex equals first and user requested close (not implemented here)
    }
    polyline = nullptr;
    actionFinished();
}

void RS_ActionDrawPolylineCustom::setRadius(double r) {
    if (r < 0.0) return;
    radius = r;
}

double RS_ActionDrawPolylineCustom::getRadius() const {
    return radius;
}

void RS_ActionDrawPolylineCustom::setVariableRadiusMode(bool on) {
    variableRadius = on;
}

bool RS_ActionDrawPolylineCustom::variableRadiusMode() const {
    return variableRadius;
}

/*
 * Geometry helper:
 *
 * Given prev point Pprev, center vertex Pcenter, and next point Pnext, compute tangent points T1 and T2
 * on the segments Pprev-Pcenter and Pcenter-Pnext respectively that are at distance from the corner
 * so an arc of radius 'radius' can fit, and compute arc center and orientation.
 *
 * Returns true if fillet can be computed and doesn't degenerate for given radius.
 */
bool RS_ActionDrawPolylineCustom::computeFilletTangents(const RS_Vector& Pprev,
                                                        const RS_Vector& Pcenter,
                                                        const RS_Vector& Pnext,
                                                        double radius,
                                                        RS_Vector& T1,
                                                        RS_Vector& T2,
                                                        RS_Vector& centerArc,
                                                        double& arcAngle,
                                                        bool& arcReversed)
{
    // Vectors from center
    RS_Vector v1 = (Pprev - Pcenter);
    RS_Vector v2 = (Pnext - Pcenter);

    double l1 = v1.magnitude();
    double l2 = v2.magnitude();

    // if either segment is zero-length, can't fillet
    if (l1 < RS_TOLERANCE || l2 < RS_TOLERANCE) {
        return false;
    }

    // unit vectors
    RS_Vector u1 = v1 / l1;
    RS_Vector u2 = v2 / l2;

    // angle between reversed u1 and u2 (internal angle at center)
    double dot = (-u1).dot(u2);
    // clamp
    dot = qBound(-1.0, dot, 1.0);
    double angle = acos(dot); // internal angle (0..pi)
    // if angle is nearly 0 or pi, fillet not meaningful
    if (angle < 1e-6 || fabs(angle - M_PI) < 1e-6) {
        return false;
    }

    // distance from corner along each segment to tangent point:
    // distance = radius / tan(angle/2)
    double dist = radius / tan(angle / 2.0);

    // if dist longer than either adjacent segment length, fillet cannot be created at full radius
    double maxDist = qMin(l1, l2);
    if (dist > maxDist - RS_TOLERANCE) {
        // too small segments; we cannot place tangents - caller should handle by skipping or clamping
        return false;
    }

    // compute tangent points
    T1 = Pcenter + u1 * (-dist); // along segment from center toward prev
    T2 = Pcenter + u2 * (-dist); // along segment from center toward next

    // compute arc center:
    // build normals to segments at tangent points pointing toward arc center,
    // intersection of the offset lines yields arc center.
    RS_Vector n1 = RS_Vector(-u1.y, u1.x); // left normal for segment 1
    RS_Vector n2 = RS_Vector(-u2.y, u2.x); // left normal for segment 2

    // The arc center lies at distance radius from each tangent along these normals.
    // two candidate centers: T1 + n1*radius or T1 - n1*radius depending on which side
    // choose the intersection of lines: (T1 + n1*s) and (T2 + n2*t) such that s = radius (or -radius)
    // We'll simply compute intersection of lines: L1 passing through T1 in direction n1, L2 passing through T2 in direction n2
    double denom = n1.x * n2.y - n1.y * n2.x;
    if (fabs(denom) < 1e-9) {
        return false; // parallel normals, degenerate
    }
    // solve: T1 + n1*s = T2 + n2*t
    RS_Vector rhs = T2 - T1;
    double s = (rhs.x * n2.y - rhs.y * n2.x) / denom;
    centerArc = T1 + n1 * s;

    // arcAngle: angle span from tangent1 to tangent2 as seen from centerArc
    RS_Vector cT1 = (T1 - centerArc);
    RS_Vector cT2 = (T2 - centerArc);
    double a1 = atan2(cT1.y, cT1.x);
    double a2 = atan2(cT2.y, cT2.x);
    arcAngle = a2 - a1;
    // normalize to -pi..pi
    while (arcAngle <= -M_PI) arcAngle += 2*M_PI;
    while (arcAngle > M_PI) arcAngle -= 2*M_PI;

    // determine direction of arc (reversed flag for RS_Arc constructor)
    arcReversed = (arcAngle < 0.0);

    return true;
}
