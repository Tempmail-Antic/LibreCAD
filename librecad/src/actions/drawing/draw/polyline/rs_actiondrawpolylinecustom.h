#ifndef RS_ACTIONDRAWPOLYLINECUSTOM_H
#define RS_ACTIONDRAWPOLYLINECUSTOM_H

#include "rs_actioninterface.h"
#include "rs_vector.h"
#include "rs_entitycontainer.h"
#include "rs_graphicview.h"

class RS_Polyline;

/**
 * RS_ActionDrawPolylineCustom
 *
 * Draw a polyline in "CUSTOM" mode: as vertices are added, automatically
 * insert fillets (arcs) of specified radius at intermediate vertices.
 *
 * Controls:
 *  - Default radius: stored in QSettings and modifiable during action.
 *  - Keyboard entry: numeric input accepted through command() (same as other actions).
 *  - Mouse wheel + modifier: adjust radius on the fly.
 *  - Preview of next segment and fillet shown in preview buffer.
 *
 * The result is a single RS_Polyline entity using bulge values for arcs
 * where possible (keeps continuity and DXF compatibility).
 */
class RS_ActionDrawPolylineCustom : public RS_ActionInterface {
    Q_OBJECT
public:
    enum Status {
        SetStartPoint = 0,
        SetNextPoint
    };

public:
    RS_ActionDrawPolylineCustom(RS_EntityContainer* container,
                                RS_GraphicView* graphicView);
    ~RS_ActionDrawPolylineCustom() override;

    void init(int status=0) override;

    QAction* createGUIAction(RS2::ActionType /*type*/, QObject* /*parent*/) override;

    void trigger() override;

    void mouseMoveEvent(QMouseEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void commandEvent(RS_CommandEvent* e) override;
    void mouseWheelEvent(QWheelEvent* e) override;

    bool updateMouseCursor() override;

    void finishPolyline();

    // set/get radius
    void setRadius(double r);
    double getRadius() const;

    // toggle variable radius mode
    void setVariableRadiusMode(bool on);
    bool variableRadiusMode() const;

    // static helper to compute fillet tangent points & bulge
    // returns true if fillet can be created for |radius|
    static bool computeFilletTangents(const RS_Vector& Pprev,
                                      const RS_Vector& Pcenter,
                                      const RS_Vector& Pnext,
                                      double radius,
                                      RS_Vector& T1,
                                      RS_Vector& T2,
                                      RS_Vector& centerArc,
                                      double& arcAngle,
                                      bool& arcReversed);

private:
    RS_Polyline* polyline;
    RS_Vector startPoint;
    RS_Vector lastPoint;
    Status status;
    double radius;            // current radius in drawing units
    bool variableRadius;      // allow per-vertex radii (toggleable)
    QString radiusTextBuffer; // typed-in radius as string (command line)
    QSettings settings;
};

#endif // RS_ACTIONDRAWPOLYLINECUSTOM_H
