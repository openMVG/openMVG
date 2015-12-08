// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <QGraphicsItem>
#include <QPointF>

class GraphWidget;
class QGraphicsSceneMouseEvent;

// Graphical movable QtGraphicItem to represent a control_point image observation
// A dynamic update of the control_point observation coordinates is performed thanks to variable reference.
class control_point_2DNode : public QGraphicsItem
{
public:
  control_point_2DNode(const QPointF& pos, double & x, double & y, size_t id_control_point);

  QRectF boundingRect() const;
  QPainterPath shape() const;
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

protected:
  QVariant itemChange(GraphicsItemChange change, const QVariant &value);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

private:
  size_t _id_control_point;
  double & _x;
  double & _y;
};

